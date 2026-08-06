from __future__ import annotations

import base64
import json
import os
import random
import sys
import time
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast, ClassVar, TypeAlias

from openai import AuthenticationError, InternalServerError, OpenAI, RateLimitError
from openai.types.chat import ChatCompletion
from openai.types.chat.chat_completion_content_part_image_param import (
    ChatCompletionContentPartImageParam,
    ImageURL,
)
from openai.types.chat.chat_completion_content_part_param import (
    ChatCompletionContentPartParam,
)
from openai.types.chat.chat_completion_content_part_text_param import (
    ChatCompletionContentPartTextParam,
)
from openai.types.chat.chat_completion_message_param import ChatCompletionMessageParam
from openai.types.chat.chat_completion_system_message_param import (
    ChatCompletionSystemMessageParam,
)
from openai.types.chat.chat_completion_user_message_param import (
    ChatCompletionUserMessageParam,
)

MessageContentItem: TypeAlias = ChatCompletionContentPartParam
ContextFile: TypeAlias = tuple[str, str]

MAX_RETRIES = 3
MAX_DELAY = 900


def resolve_api_key(server_type: str) -> str | None:
    """Resolve the API key based on server type."""
    if server_type == "openai":
        key = os.getenv("OPENAI_API_KEY")
        if not key:
            raise ValueError("OpenAI API key is not provided in credentials!")
        return key
    elif server_type == "custom":
        key = os.getenv("CUSTOM_SERVER_API_KEY")
        return key if key else None
    else:
        raise ValueError(f"Unknown server type: {server_type}")


def resolve_base_url(server_type: str) -> str | None:
    """Resolve the base URL based on server type."""
    if server_type == "custom":
        url = os.getenv("CUSTOM_SERVER_URL")
        if not url:
            raise ValueError("Custom server URL is not provided in credentials!")
        if not url.startswith(("http://", "https://")):
            raise ValueError(
                f"Custom server URL must start with http:// or https://, got: {url}"
            )
        return url
    return None


def build_client(base_url: str | None, api_key: str | None) -> OpenAI:
    """Create an OpenAI client, optionally with a custom base URL."""
    kwargs: dict = {}
    if api_key:
        kwargs["api_key"] = api_key
    elif base_url:
        # Custom servers (e.g. Ollama, vLLM) may not require authentication.
        # The OpenAI SDK requires an api_key, so use a placeholder.
        kwargs["api_key"] = "not-needed"
    if base_url:
        kwargs["base_url"] = base_url
    return OpenAI(**kwargs)


@dataclass(frozen=True)
class MessageBuilder:
    question: str
    context_files: Sequence[ContextFile]

    _MEDIA_TYPE_MAP: ClassVar[dict[str, str]] = {
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".png": "image/png",
        ".gif": "image/gif",
        ".webp": "image/webp",
    }

    _MAX_IMAGE_BYTES: ClassVar[int] = 20 * 1024 * 1024

    def build(self) -> list[MessageContentItem]:
        """Construct the completion request payload."""
        message: list[MessageContentItem] = [{"type": "text", "text": self.question}]

        for path, file_type in self.context_files:
            if file_type == "image":
                message.append(self._build_image_content(path))
            else:
                message.append(self._build_text_content(path))

        return message

    def _build_image_content(self, path: str) -> ChatCompletionContentPartImageParam:
        """Encode an image context file for model consumption."""
        if os.path.getsize(path) > self._MAX_IMAGE_BYTES:
            raise ValueError(
                f"File {path} exceeds the 20MB limit and will not be processed."
            )

        _, ext = os.path.splitext(path)
        media_type = self._MEDIA_TYPE_MAP.get(ext.lower(), "image/jpeg")
        with open(path, "rb") as img_file:
            image_data = base64.standard_b64encode(img_file.read()).decode("utf-8")

        image_url_payload = ImageURL(
            url=f"data:{media_type};base64,{image_data}", detail="auto"
        )
        return ChatCompletionContentPartImageParam(
            type="image_url", image_url=image_url_payload
        )

    def _build_text_content(self, path: str) -> ChatCompletionContentPartTextParam:
        """Read a text context file and wrap it in a templated message."""
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as text_file:
                file_content = text_file.read()
        except OSError as exc:
            raise ValueError(f"Error reading file {path}: {exc}") from exc

        basename = os.path.basename(path)
        return ChatCompletionContentPartTextParam(
            type="text",
            text=f"--- Content of {basename} ---\n{file_content}\n",
        )


def parse_context_files(raw: str) -> list[ContextFile]:
    """Parse and validate the JSON encoded context file descriptors."""
    try:
        decoded = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError("Invalid JSON payload for context files.") from exc

    if not isinstance(decoded, list):
        raise ValueError("Context files payload must be a list.")

    parsed: list[ContextFile] = []
    for entry in decoded:
        if (
            isinstance(entry, list)
            and len(entry) == 2
            and isinstance(entry[0], str)
            and isinstance(entry[1], str)
        ):
            parsed.append((entry[0], entry[1]))
        else:
            raise ValueError(
                "Each context file entry must be a pair of strings [path, type]."
            )

    return parsed


def build_messages(
    question: str,
    context_files: Sequence[ContextFile],
    system_message: str | None = None,
) -> list[ChatCompletionMessageParam]:
    """Build the full message list including optional system message and user content."""
    message_content = MessageBuilder(
        question=question, context_files=context_files
    ).build()

    messages: list[ChatCompletionMessageParam] = []

    if system_message:
        messages.append(
            cast(
                ChatCompletionMessageParam,
                ChatCompletionSystemMessageParam(
                    role="system", content=system_message
                ),
            )
        )

    user_message = ChatCompletionUserMessageParam(role="user", content=list(message_content))
    messages.append(cast(ChatCompletionMessageParam, user_message))
    return messages


def call_chat_completion(
    client: OpenAI,
    model: str,
    messages: list[ChatCompletionMessageParam],
    temperature: float | None = None,
    max_tokens: int | None = None,
    top_p: float | None = None,
) -> ChatCompletion:
    """Request a chat completion with optional parameters."""
    api_params: dict = {"model": model, "messages": messages}
    if temperature is not None:
        api_params["temperature"] = temperature
    if max_tokens is not None:
        api_params["max_tokens"] = max_tokens
    if top_p is not None:
        api_params["top_p"] = top_p

    return client.chat.completions.create(**api_params)


def _call_with_retries(
    client: OpenAI,
    model: str,
    messages: list[ChatCompletionMessageParam],
    temperature: float | None,
    max_tokens: int | None,
    top_p: float | None,
) -> ChatCompletion | None:
    """Call chat completion with exponential backoff retry on server errors."""
    for attempt in range(MAX_RETRIES):
        try:
            return call_chat_completion(
                client, model, messages, temperature, max_tokens, top_p
            )
        except (InternalServerError, RateLimitError) as exc:
            if attempt == MAX_RETRIES - 1:
                print(f"Max retries reached. Last error: {exc}")
                return None
            sleep_time = min(2**attempt + random.uniform(0, 1), MAX_DELAY)
            print(
                f"Server error encountered ({exc}). Retrying in {sleep_time:.2f}s..."
            )
            time.sleep(sleep_time)
        except AuthenticationError as exc:
            print(f"Authentication error: {exc}")
            return None
        except Exception as exc:  # noqa: BLE001 - keep reporting unexpected errors
            print(f"An error occurred: {exc}")
            return None
    return None


def main(argv: Sequence[str]) -> int:
    if len(argv) < 9:
        print(
            "Usage: chatgpt.py <context_files_json> <question> <model> "
            "<server_type> <temperature> <max_tokens> <top_p> <system_message>"
        )
        return 1

    try:
        context_files = parse_context_files(argv[1])
    except ValueError as exc:
        print(str(exc))
        return 1

    question = argv[2].replace("__cn__", "\n")
    model = argv[3]
    server_type = argv[4]
    temperature_arg = argv[5]
    max_tokens_arg = argv[6]
    top_p_arg = argv[7]
    system_message_arg = argv[8]

    temperature = float(temperature_arg) if temperature_arg and temperature_arg != "None" else None
    max_tokens = int(max_tokens_arg) if max_tokens_arg and max_tokens_arg != "None" else None
    top_p = float(top_p_arg) if top_p_arg and top_p_arg != "None" else None
    system_message = system_message_arg.replace("__cn__", "\n") if system_message_arg and system_message_arg != "None" else None

    try:
        api_key = resolve_api_key(server_type)
    except ValueError as exc:
        print(str(exc))
        return 1

    try:
        base_url = resolve_base_url(server_type)
    except ValueError as exc:
        print(str(exc))
        return 1

    try:
        client = build_client(base_url, api_key)
    except Exception as exc:  # noqa: BLE001
        print(f"An error occurred: {exc}")
        return 1

    try:
        messages = build_messages(question, context_files, system_message)
    except ValueError as exc:
        print(str(exc))
        return 1

    response = _call_with_retries(client, model, messages, temperature, max_tokens, top_p)
    if response is None:
        return 1

    if not response.choices:
        server_label = server_type if server_type != "openai" else "OpenAI"
        print(
            f"No output was generated!\n"
            f"Please ensure that your {server_label} account has sufficient credits "
            f"or that the model '{model}' is available on the configured server."
        )
        return 1

    message = response.choices[0].message
    content = getattr(message, "content", None)
    if not content:
        server_label = server_type if server_type != "openai" else "OpenAI"
        print(
            f"No output was generated!\n"
            f"Please ensure that your {server_label} account has sufficient credits "
            f"or that the model '{model}' is available on the configured server."
        )
        return 1

    print(
        f"Successfully generated response for:\n{question[:100]}{'...' if len(question) > 100 else ''}"
    )
    with open("output.md", "w", encoding="utf-8") as file_handle:
        file_handle.write(content)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
