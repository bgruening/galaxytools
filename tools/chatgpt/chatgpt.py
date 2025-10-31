from __future__ import annotations

import base64
import json
import os
import sys
from dataclasses import dataclass
from typing import ClassVar, Iterable, List, Sequence, Tuple, cast

from openai import AuthenticationError, OpenAI
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
from openai.types.chat.chat_completion_user_message_param import (
    ChatCompletionUserMessageParam,
)


MessageContentItem = ChatCompletionContentPartParam
ContextFile = Tuple[str, str]


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

    def build(self) -> List[MessageContentItem]:
        """Construct the completion request payload."""
        message: List[MessageContentItem] = [{"type": "text", "text": self.question}]

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


def parse_context_files(raw: str) -> List[ContextFile]:
    """Parse and validate the JSON encoded context file descriptors."""
    try:
        decoded = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError("Invalid JSON payload for context files.") from exc

    if not isinstance(decoded, list):
        raise ValueError("Context files payload must be a list.")

    parsed: List[ContextFile] = []
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
    question: str, context_files: Sequence[ContextFile]
) -> List[MessageContentItem]:
    """Helper to hide the dataclass implementation detail from callers."""
    return MessageBuilder(question=question, context_files=context_files).build()


def call_chat_completion(
    client: OpenAI, model: str, messages: Iterable[MessageContentItem]
) -> ChatCompletion:
    """Request a chat completion using the given client."""
    user_message = ChatCompletionUserMessageParam(role="user", content=list(messages))
    payload: List[ChatCompletionMessageParam] = [
        cast(ChatCompletionMessageParam, user_message)
    ]
    # Some models (e.g. the GPT-4.1 family) reject non-default temperature values,
    # so rely on the API default rather than forcing a custom value here.
    return client.chat.completions.create(
        model=model,
        messages=payload,
    )


def main(argv: Sequence[str]) -> int:
    if len(argv) < 4:
        print("Usage: chatgpt.py <context_files_json> <question> <model>")
        return 1

    try:
        context_files = parse_context_files(argv[1])
    except ValueError as exc:
        print(str(exc))
        return 1

    question = argv[2]
    model = argv[3]

    openai_api_key = os.getenv("OPENAI_API_KEY")
    if not openai_api_key:
        print("OpenAI API key is not provided in credentials!")
        return 1

    client = OpenAI(api_key=openai_api_key)

    try:
        message_content = build_messages(question, context_files)
        response = call_chat_completion(client, model, message_content)
    except AuthenticationError as exc:
        print(f"Authentication error: {exc}")
        return 1
    except ValueError as exc:
        print(str(exc))
        return 1
    except Exception as exc:  # noqa: BLE001 - keep reporting unexpected OpenAI errors
        print(f"An error occurred: {exc}")
        return 1

    if not response.choices:
        print(
            "No output was generated!\nPlease ensure that your OpenAI account has sufficient credits.\n"
            "You can check your balance here: https://platform.openai.com/settings/organization/billing"
        )
        return 1

    message = response.choices[0].message
    content = getattr(message, "content", None)
    if not content:
        print(
            "No output was generated!\nPlease ensure that your OpenAI account has sufficient credits.\n"
            "You can check your balance here: https://platform.openai.com/settings/organization/billing"
        )
        return 1

    print("Output has been saved!")
    with open("output.txt", "w", encoding="utf-8") as file_handle:
        file_handle.write(content)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
