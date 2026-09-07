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

from openai import (
    APIConnectionError,
    AuthenticationError,
    InternalServerError,
    OpenAI,
    RateLimitError,
)
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

BILLING_URL = "https://platform.openai.com/settings/organization/billing"

MAX_ERROR_CHARS = 300

REQUEST_TIMEOUT = 300.0


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
                "Custom server URL must start with http:// or https://"
            )
        # The SDK appends "chat/completions" to the base URL's path *and* its
        # query string, so a trailing "?" makes the query swallow the suffix
        # and leaves the request path entirely under the credential's control
        # -- turning an LLM endpoint setting into an arbitrary POST target.
        if "?" in url or "#" in url:
            raise ValueError(
                "Custom server URL must not contain a query string or fragment."
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
    # Retry once, here, rather than letting the SDK multiply each attempt.
    # The timeout is explicit because the SDK's 600s default, multiplied by
    # this tool's retries, can hold a Galaxy job slot for half an hour.
    return OpenAI(max_retries=0, timeout=REQUEST_TIMEOUT, **kwargs)


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
        try:
            size = os.path.getsize(path)
        except OSError as exc:
            raise ValueError(f"Error reading file {path}: {exc}") from exc

        if size > self._MAX_IMAGE_BYTES:
            raise ValueError(
                f"File {path} exceeds the 20MB limit and will not be processed."
            )

        _, ext = os.path.splitext(path)
        media_type = self._MEDIA_TYPE_MAP.get(ext.lower(), "image/jpeg")
        try:
            with open(path, "rb") as img_file:
                image_data = base64.standard_b64encode(img_file.read()).decode("utf-8")
        except OSError as exc:
            raise ValueError(f"Error reading file {path}: {exc}") from exc

        image_url_payload = ImageURL(
            url=f"data:{media_type};base64,{image_data}", detail="auto"
        )
        return ChatCompletionContentPartImageParam(
            type="image_url", image_url=image_url_payload
        )

    def _build_text_content(self, path: str) -> ChatCompletionContentPartTextParam:
        """Read a text context file and wrap it in a templated message."""
        file_content = read_text_file(path)

        basename = os.path.basename(path)
        return ChatCompletionContentPartTextParam(
            type="text",
            text=f"--- Content of {basename} ---\n{file_content}\n",
        )


def read_text_file(path: str) -> str:
    """Read a UTF-8 text file, reporting a tool-level error on failure."""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as text_file:
            return text_file.read()
    except OSError as exc:
        raise ValueError(f"Error reading file {path}: {exc}") from exc


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


# OpenAI models offered by this tool that accept temperature/top_p on Chat
# Completions. Determined empirically against the live API on 2026-09-07: every
# other model in the option list answers a request carrying temperature with
# 400 unsupported_value, and one carrying top_p with 400 unsupported_parameter.
# gpt-5.4 accepts both because its default reasoning effort is "none".
# A prefix rule cannot express this -- "gpt-6-astra" does not start with
# "gpt-5", and "gpt-5.4" does -- so keep an explicit set, and re-check it when
# adding an option. Anything not listed is treated as refusing the parameters:
# sending one to a model that refuses it fails the job, omitting it only
# costs a note.
SAMPLING_MODELS = frozenset({"gpt-4.1", "gpt-4o", "gpt-5.4"})


def uses_fixed_sampling(model: str) -> bool:
    """Whether an OpenAI model refuses temperature/top_p."""
    return model not in SAMPLING_MODELS


def build_api_params(
    model: str,
    messages: list[ChatCompletionMessageParam],
    server_type: str,
    temperature: float | None = None,
    max_tokens: int | None = None,
    top_p: float | None = None,
) -> dict:
    """Assemble the request parameters, adapted to the target server."""
    api_params: dict = {"model": model, "messages": messages}
    fixed_sampling = server_type == "openai" and uses_fixed_sampling(model)

    for name, value in (("temperature", temperature), ("top_p", top_p)):
        if value is None:
            continue
        if fixed_sampling:
            if value != 1.0:
                print(
                    f"Note: {name} is not sent for '{model}' -- reasoning "
                    f"models on the OpenAI API reject it; the requested value "
                    f"{value} was ignored."
                )
            continue
        api_params[name] = value

    if max_tokens is not None:
        # OpenAI deprecated ``max_tokens`` and rejects it outright on the gpt-5
        # family; ``max_completion_tokens`` is accepted by every current OpenAI
        # model. Custom servers are the mirror image -- Ollama's OpenAI shim
        # still only understands ``max_tokens``.
        if server_type == "openai":
            api_params["max_completion_tokens"] = max_tokens
        else:
            api_params["max_tokens"] = max_tokens

    return api_params


def call_chat_completion(client: OpenAI, api_params: dict) -> ChatCompletion:
    """Request a chat completion with the prepared parameters."""
    return client.chat.completions.create(**api_params)


def describe_error(exc: Exception, trusted: bool = False) -> str:
    """Summarise an API error without echoing a server's response back.

    The custom server URL comes from the user, so the job node can be pointed
    at any host it can reach and its reply must not reach the job log. Only
    api.openai.com is ``trusted``: for anything else the SDK may hand back the
    whole response body -- it only unwraps an "error" key when one is present
    -- so the body could be any internal service's, and none of it is shown.
    """
    parts: list[str] = []

    status = getattr(exc, "status_code", None)
    if status:
        parts.append(f"HTTP {status}")

    body = getattr(exc, "body", None)
    if isinstance(body, dict):
        if trusted:
            code = body.get("code")
            if isinstance(code, (str, int)) and str(code).strip():
                parts.append(f"code {str(code).strip()[:MAX_ERROR_CHARS]}")
            message = body.get("message")
            if isinstance(message, str) and message.strip():
                parts.append(message.strip()[:MAX_ERROR_CHARS])
        else:
            parts.append("the server returned an error response")
    elif body is not None:
        parts.append("the server returned an unexpected non-JSON response")

    if not parts:
        # No HTTP response at all. The cause is raised locally by httpx or the
        # OS, so it is safe to show, and it is the only thing that separates a
        # refused connection from a DNS, TLS or timeout failure -- the most
        # likely outcome of a mistyped custom server URL.
        detail = ""
        if isinstance(exc, APIConnectionError) and exc.__cause__ is not None:
            cause = exc.__cause__
            detail = f" ({type(cause).__name__}: {cause})"[:MAX_ERROR_CHARS]
        parts.append(f"{type(exc).__name__}{detail}")
    return "; ".join(parts)


def is_quota_error(exc: Exception) -> bool:
    """Whether a rate limit error is really an exhausted account balance."""
    for attr in ("code", "type"):
        if getattr(exc, attr, None) == "insufficient_quota":
            return True
    if isinstance(getattr(exc, "body", None), dict):
        return False
    return "insufficient_quota" in str(exc)


def _call_with_retries(
    client: OpenAI,
    api_params: dict,
    server_type: str,
) -> ChatCompletion | None:
    """Call chat completion with exponential backoff retry on server errors."""
    for attempt in range(MAX_RETRIES):
        try:
            return call_chat_completion(client, api_params)
        except (APIConnectionError, InternalServerError, RateLimitError) as exc:
            # An exhausted balance is reported as a rate limit, but retrying it
            # is pointless and hides a billing problem behind a server error.
            trusted = server_type == "openai"
            if is_quota_error(exc):
                if server_type == "openai":
                    print(
                        "Insufficient quota!\n"
                        "Please ensure that your OpenAI account has sufficient credits.\n"
                        f"You can check your balance here: {BILLING_URL}"
                    )
                else:
                    print(
                        "Insufficient quota reported by the configured server: "
                        f"{describe_error(exc, trusted)}"
                    )
                return None
            if attempt == MAX_RETRIES - 1:
                print(f"Max retries reached. Last error: {describe_error(exc, trusted)}")
                return None
            sleep_time = 2**attempt + random.uniform(0, 1)
            print(
                f"Server error encountered ({describe_error(exc, trusted)}). "
                f"Retrying in {sleep_time:.2f}s..."
            )
            time.sleep(sleep_time)
        except AuthenticationError as exc:
            print(f"Authentication error: {describe_error(exc, server_type == 'openai')}")
            return None
        except Exception as exc:  # noqa: BLE001 - keep reporting unexpected errors
            print(f"An error occurred: {describe_error(exc, server_type == 'openai')}")
            return None
    return None


def main(argv: Sequence[str]) -> int:
    if len(argv) < 9:
        print(
            "Usage: chatgpt.py <context_files_json> <prompt_file> <model> "
            "<server_type> <temperature> <max_tokens> <top_p> <system_message_file>"
        )
        return 1

    try:
        context_files = parse_context_files(argv[1])
    except ValueError as exc:
        print(str(exc))
        return 1

    model = argv[3]
    server_type = argv[4]
    temperature_arg = argv[5]
    max_tokens_arg = argv[6]
    top_p_arg = argv[7]

    # The prompt and the system message are passed as files rather than on the
    # command line so that Galaxy's parameter sanitizer can be turned off for
    # them: on the command line an apostrophe would break the shell quoting and
    # every non-ASCII character would be replaced with a literal "X".
    try:
        question = read_text_file(argv[2])
        system_message = read_text_file(argv[8]).strip() or None
    except ValueError as exc:
        print(str(exc))
        return 1

    if not question.strip():
        print("The prompt is empty!")
        return 1

    temperature = float(temperature_arg) if temperature_arg and temperature_arg != "None" else None
    max_tokens = int(max_tokens_arg) if max_tokens_arg and max_tokens_arg != "None" else None
    top_p = float(top_p_arg) if top_p_arg and top_p_arg != "None" else None

    try:
        api_key = resolve_api_key(server_type)
        base_url = resolve_base_url(server_type)
    except ValueError as exc:
        print(str(exc))
        return 1

    try:
        client = build_client(base_url, api_key)
    except Exception:  # noqa: BLE001
        print(
            "The configured server URL could not be used to build a client; "
            "check its host and port."
        )
        return 1

    try:
        messages = build_messages(question, context_files, system_message)
    except ValueError as exc:
        print(str(exc))
        return 1

    api_params = build_api_params(
        model, messages, server_type, temperature, max_tokens, top_p
    )
    response = _call_with_retries(client, api_params, server_type)
    if response is None:
        return 1

    choice = response.choices[0] if response.choices else None
    message = getattr(choice, "message", None)
    content = getattr(message, "content", None)
    if content is not None and not isinstance(content, str):
        print(
            "The server returned a response in an unexpected format; this tool "
            "expects an OpenAI-compatible chat completion."
        )
        return 1
    if not content:
        refusal = getattr(message, "refusal", None)
        if refusal:
            print(f"The model declined to answer: {str(refusal)[:MAX_ERROR_CHARS]}")
        elif choice is not None and choice.finish_reason == "length":
            print(
                "No output was generated!\n"
                "The response hit the 'Max tokens' limit before any answer was "
                "produced. On the gpt-5 models that budget also covers hidden "
                "reasoning tokens, so raise Max tokens or leave it unset."
            )
        elif server_type == "openai":
            print(
                "No output was generated!\n"
                "Please ensure that your OpenAI account has sufficient credits "
                f"or that the model '{model}' is available.\n"
                f"You can check your balance here: {BILLING_URL}"
            )
        else:
            print(
                "No output was generated!\n"
                f"Please ensure that the model '{model}' is available on the "
                "configured server."
            )
        return 1

    with open("output.md", "w", encoding="utf-8") as file_handle:
        file_handle.write(content)

    print(
        f"Successfully generated response for:\n{question[:100]}{'...' if len(question) > 100 else ''}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
