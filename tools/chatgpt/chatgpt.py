from __future__ import annotations

import base64
import json
import os
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from typing import ClassVar, TypeAlias

from openai import (
    APIConnectionError,
    APIStatusError,
    APITimeoutError,
    DefaultHttpxClient,
    OpenAI,
)
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

MessageContentItem: TypeAlias = ChatCompletionContentPartParam
ContextFile: TypeAlias = tuple[str, str]

BILLING_URL = "https://platform.openai.com/settings/organization/billing"
MAX_ERROR_CHARS = 300

# Options a model may refuse, e.g. reasoning models refuse temperature.
OPTIONAL_PARAMS = ("temperature", "top_p")

# Servers word their errors differently, so look for these words in the error.
REFUSED = re.compile(r"unsupported|not supported|not support|n't support|extra|unrecognized|must be|not permitted")
TOO_LONG = re.compile(r"context.?(length|window|size)|maximum context|too large for model|longer than the model|reduce the length|max_new_tokens")
NO_IMAGES = re.compile(r"(image|multimodal)[^.]{0,40}(not (supported|enabled)|only supported)|not a multimodal|support image input|does not support images")
NO_MODEL = re.compile(r"model[^.]{0,60}(not found|does not exist|not available)|invalid model|not a valid model|model_not_found|deploymentnotfound")
NO_CREDIT = re.compile(r"quota|budget|credit|billing")


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


def make_client(server_type: str) -> OpenAI:
    """Create the client from the credentials Galaxy puts in the environment."""
    if server_type == "openai":
        api_key = os.getenv("OPENAI_API_KEY", "").strip()
        if not api_key:
            raise ValueError("OpenAI API key is not provided in credentials!")
        base_url = None
    else:
        base_url = os.getenv("CUSTOM_SERVER_URL", "").strip()
        if not base_url:
            raise ValueError("Custom server URL is not provided in credentials!")
        # The SDK adds "chat/completions" after the URL's path and query, so a
        # "?" or "#" would let the URL send the request to any path.
        if not base_url.lower().startswith(("http://", "https://")) or "?" in base_url or "#" in base_url:
            raise ValueError("Error: The server URL must start with http:// or https:// and must not contain '?' or '#'.")
        # Many local servers need no key, but the SDK needs one.
        api_key = os.getenv("CUSTOM_SERVER_API_KEY", "").strip() or "not-needed"
    # No redirects: they could send the request, key included, to another URL.
    return OpenAI(api_key=api_key, base_url=base_url, http_client=DefaultHttpxClient(follow_redirects=False))


@dataclass
class Reply:
    content: str
    refusal: str
    finish_reason: str | None


def call_chat_completion(client: OpenAI, api_params: dict) -> Reply | None:
    """Request a chat completion and collect the streamed answer.

    Streaming keeps the connection busy while a reasoning model thinks, so a
    long answer does not hit the read timeout. Returns None if the reply is
    not a chat answer at all, e.g. a web page from a wrong URL.
    """
    content, refusal, finish_reason, chunks = [], [], None, 0
    with client.chat.completions.create(**api_params, stream=True) as stream:
        for chunk in stream:
            chunks += 1
            if not getattr(chunk, "choices", None):
                continue
            choice = chunk.choices[0]
            if isinstance(choice.delta.content, str):
                content.append(choice.delta.content)
            if isinstance(getattr(choice.delta, "refusal", None), str):
                refusal.append(choice.delta.refusal)
            finish_reason = choice.finish_reason or finish_reason
    if not chunks:
        return None
    return Reply("".join(content), "".join(refusal), finish_reason)


def error_text(exc: Exception) -> str:
    """All text of an error, in lower case, to look for the words above."""
    body = getattr(exc, "body", None)
    body_text = json.dumps(body) if isinstance(body, (dict, list)) else str(body or "")
    return f"{getattr(exc, 'param', None) or ''} {exc} {body_text}".lower()


def refused_params(exc: APIStatusError, api_params: dict) -> list[str]:
    """The optional parameters the server refused, if that is the error."""
    text = error_text(exc)
    if exc.status_code not in (400, 422) or TOO_LONG.search(text):
        return []
    sent = [p for p in OPTIONAL_PARAMS if p in api_params]
    if getattr(exc, "param", None) in sent:
        return [exc.param]
    if exc.status_code == 422 or REFUSED.search(text):
        return [p for p in sent if p in text]
    return []


def explain_error(exc: Exception, model: str, server_type: str) -> str:
    """A short, plain message for a failed request.

    The custom server URL comes from the user, so the job node can be pointed
    at any host it can reach. Its reply is only searched for known words and
    never written to the job log; only OpenAI's own error message is shown.
    """
    status = getattr(exc, "status_code", None)
    text = error_text(exc)
    if isinstance(exc, APITimeoutError):
        message = "The server took too long to answer. Try again later."
    elif isinstance(exc, APIConnectionError):
        message = "Could not reach the server. Check the server URL."
    elif status == 401:
        message = "The API key was not accepted. Check the key in your credentials."
    elif status == 403:
        message = "Access denied. Your key may not be allowed to use this model."
    elif TOO_LONG.search(text):
        message = "The prompt and context files are too long for this model. Use fewer or smaller files."
    elif NO_IMAGES.search(text):
        message = "This model cannot read images. Choose a model that can, or remove the images."
    elif NO_MODEL.search(text) or (status == 404 and "model" in text):
        message = f"The server does not know the model '{model}'. Check the model name."
    elif status in (402, 429) and NO_CREDIT.search(text):
        message = "Your account has no credits or budget left."
        if server_type == "openai":
            message += f" Check your balance: {BILLING_URL}"
    elif status == 429:
        message = "Too many requests. Wait a few minutes and try again."
    elif status and status >= 500:
        message = "The server had a problem. Try again later."
    elif status in (404, 405) or (status and 300 <= status < 400):
        message = "The server URL seems wrong. Check that it ends with the API path, for example /v1 or /api."
    else:
        message = "The server could not handle the request."

    # Without a status, show only the error type: a connection error's
    # message can contain the request headers, including the API key.
    details = f"HTTP {status}" if status else type(exc.__cause__ or exc).__name__
    body = getattr(exc, "body", None)
    if server_type == "openai" and isinstance(body, dict) and isinstance(body.get("message"), str):
        details += f" - {body['message'][:MAX_ERROR_CHARS]}"
    return f"Error: {message}\nDetails: {details}"


def request_completion(client: OpenAI, api_params: dict, server_type: str) -> Reply | None:
    """Call chat completion and report any error. The SDK handles retries.

    If the model refuses temperature or top_p, send the request again
    without them and say so in the job log.
    """
    model = api_params["model"]
    try:
        reply = call_chat_completion(client, api_params)
        if reply is None:
            print(
                "Error: The server URL seems wrong. Its reply is not a chat answer. "
                "Check that it ends with the API path, for example /v1 or /api."
            )
        return reply
    except APIStatusError as exc:
        refused = refused_params(exc, api_params)
        if refused:
            for param in refused:
                print(f"Note: model '{model}' does not accept {param}, so it was left out.")
            api_params = {k: v for k, v in api_params.items() if k not in refused}
            return request_completion(client, api_params, server_type)
        print(explain_error(exc, model, server_type))
    except Exception as exc:  # noqa: BLE001 - report every failure plainly
        print(explain_error(exc, model, server_type))
    return None


def optional(value: str, kind: type) -> float | int | None:
    """Galaxy passes an unset optional number as an empty string or 'None'."""
    return kind(value) if value not in ("", "None") else None


def main(argv: Sequence[str]) -> int:
    if len(argv) < 9:
        print(
            "Usage: chatgpt.py <context_files_json> <prompt_file> <model> "
            "<server_type> <temperature> <max_tokens> <top_p> <system_message_file>"
        )
        return 1

    model, server_type = argv[3], argv[4]
    # The prompt and the system message come as files, so Galaxy does not
    # need to sanitize them (that would change quotes and non-ASCII text).
    with open(argv[2], encoding="utf-8") as file_handle:
        question = file_handle.read()
    with open(argv[8], encoding="utf-8") as file_handle:
        system_message = file_handle.read().strip()
    if not question.strip():
        print("The prompt is empty!")
        return 1

    try:
        client = make_client(server_type)
        messages: list[ChatCompletionMessageParam] = []
        if system_message:
            messages.append({"role": "system", "content": system_message})
        content = MessageBuilder(question, parse_context_files(argv[1])).build()
        messages.append({"role": "user", "content": content})
    except ValueError as exc:
        print(str(exc))
        return 1
    except Exception:  # noqa: BLE001 - e.g. a URL the SDK cannot parse
        print("Error: The server URL is not valid. Check it in your credentials.")
        return 1

    api_params: dict = {"model": model, "messages": messages}
    for name, value in (("temperature", optional(argv[5], float)), ("top_p", optional(argv[7], float))):
        if value is not None:
            api_params[name] = value
    max_tokens = optional(argv[6], int)
    if max_tokens is not None:
        # OpenAI's reasoning models only take max_completion_tokens, while
        # Ollama only takes max_tokens.
        api_params["max_completion_tokens" if server_type == "openai" else "max_tokens"] = max_tokens

    reply = request_completion(client, api_params, server_type)
    if reply is None:
        return 1
    if not reply.content:
        if reply.refusal:
            print(f"The model declined to answer: {reply.refusal[:MAX_ERROR_CHARS]}")
        elif reply.finish_reason == "length":
            print(
                "Error: The answer was empty because it hit the Max tokens limit. "
                "Raise Max tokens or leave it empty. Reasoning models also "
                "count their hidden thinking in this limit."
            )
        else:
            print(f"Error: The model '{model}' gave an empty answer. Try again or choose another model.")
        return 1

    with open("output.md", "w", encoding="utf-8") as file_handle:
        file_handle.write(reply.content)
    if reply.finish_reason in ("length", "content_filter"):
        print(f"Warning: the answer may be incomplete (finish_reason: {reply.finish_reason}).")
    print(
        f"Successfully generated response for:\n{question[:100]}{'...' if len(question) > 100 else ''}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
