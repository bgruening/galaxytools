from __future__ import annotations

import base64
import json
import os
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast, ClassVar, TypeAlias

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
from openai.types.chat.chat_completion_system_message_param import (
    ChatCompletionSystemMessageParam,
)
from openai.types.chat.chat_completion_user_message_param import (
    ChatCompletionUserMessageParam,
)

MessageContentItem: TypeAlias = ChatCompletionContentPartParam
ContextFile: TypeAlias = tuple[str, str]

BILLING_URL = "https://platform.openai.com/settings/organization/billing"

MAX_ERROR_CHARS = 300


def resolve_api_key(server_type: str) -> str | None:
    """Resolve the API key based on server type."""
    if server_type == "openai":
        key = os.getenv("OPENAI_API_KEY", "").strip()
        if not key:
            raise ValueError("OpenAI API key is not provided in credentials!")
        return key
    elif server_type == "custom":
        key = os.getenv("CUSTOM_SERVER_API_KEY", "").strip()
        return key or None
    else:
        raise ValueError(f"Unknown server type: {server_type}")


def resolve_base_url(server_type: str) -> str | None:
    """Resolve the base URL based on server type."""
    if server_type == "custom":
        url = os.getenv("CUSTOM_SERVER_URL", "").strip()
        if not url:
            raise ValueError("Custom server URL is not provided in credentials!")
        if not url.lower().startswith(("http://", "https://")):
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
    # Keep the SDK's default timeouts and retries, but do not follow redirects:
    # a redirect would let the server send the request, prompt and key
    # included, to any other URL.
    return OpenAI(http_client=DefaultHttpxClient(follow_redirects=False), **kwargs)


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
    if temperature is not None:
        api_params["temperature"] = temperature
    if top_p is not None:
        api_params["top_p"] = top_p

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


# Options a model may refuse, e.g. reasoning models refuse temperature.
OPTIONAL_PARAMS = ("temperature", "top_p")

# Servers word their errors differently, so look for these words in the error.
REFUSED = re.compile(r"unsupported|not supported|not support|n't support|extra|unrecognized|must be|not permitted")
TOO_LONG = re.compile(r"context.?(length|window|size)|maximum context|too large for model|longer than the model|reduce the length|max_new_tokens")
NO_IMAGES = re.compile(r"(image|multimodal)[^.]{0,40}(not (supported|enabled)|only supported)|not a multimodal|support image input|does not support images")
NO_MODEL = re.compile(r"model[^.]{0,60}(not found|does not exist|not available)|invalid model|not a valid model|model_not_found|deploymentnotfound")
NO_CREDIT = re.compile(r"quota|budget|credit|billing")


def error_text(exc: Exception) -> str:
    """All text of an error, in lower case, to look for the words above."""
    body = getattr(exc, "body", None)
    parts = [str(getattr(exc, "param", None) or ""), str(exc)]
    if isinstance(body, (dict, list)):
        parts.append(json.dumps(body))
    elif body:
        parts.append(str(body))
    return " ".join(parts).lower()


def refused_params(exc: APIStatusError, api_params: dict) -> list[str]:
    """The optional parameters the server refused, if that is the error."""
    if exc.status_code not in (400, 422):
        return []
    text = error_text(exc)
    if TOO_LONG.search(text):
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

    if status:
        details = f"HTTP {status}"
    else:
        # Only the type: the message of a connection error can contain the
        # request headers, including the API key.
        details = type(exc.__cause__ or exc).__name__
    body = getattr(exc, "body", None)
    if server_type == "openai" and isinstance(body, dict) and isinstance(body.get("message"), str):
        details += f" - {body['message'][:MAX_ERROR_CHARS]}"
    return f"Error: {message}\nDetails: {details}"


def request_completion(
    client: OpenAI,
    api_params: dict,
    server_type: str,
) -> Reply | None:
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
        print("Error: The server URL is not valid. Check it in your credentials.")
        return 1

    try:
        messages = build_messages(question, context_files, system_message)
    except ValueError as exc:
        print(str(exc))
        return 1

    api_params = build_api_params(
        model, messages, server_type, temperature, max_tokens, top_p
    )
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
            print(
                f"Error: The model '{model}' gave an empty answer. "
                "Try again or choose another model."
            )
        return 1

    with open("output.md", "w", encoding="utf-8") as file_handle:
        file_handle.write(reply.content)

    if reply.finish_reason in ("length", "content_filter"):
        print(
            "Warning: the answer may be incomplete "
            f"(finish_reason: {reply.finish_reason})."
        )
    print(
        f"Successfully generated response for:\n{question[:100]}{'...' if len(question) > 100 else ''}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
