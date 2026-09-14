import hashlib
import json
import os
import random
import sys
import time

import yaml
from openai import (
    APIConnectionError,
    APIError,
    APITimeoutError,
    BadRequestError,
    InternalServerError,
    OpenAI,
    RateLimitError,
)

context_files = json.loads(sys.argv[1])
question = sys.argv[2]
model = sys.argv[3]
model_type = sys.argv[4]
temperature_arg = sys.argv[5]
temperature = float(temperature_arg) if temperature_arg else None
provider = sys.argv[6]
# Galaxy user id, for request attribution on the proxy. The literal
# "Anonymous" for anonymous sessions, in which case no per-user id is sent.
galaxy_user_id = sys.argv[7] if len(sys.argv) > 7 else ""
if galaxy_user_id == "Anonymous":
    galaxy_user_id = ""
# Galaxy instance URL ($__galaxy_url__), used to namespace the user id below.
galaxy_url = sys.argv[8] if len(sys.argv) > 8 else ""

litellm_config_file = os.environ.get("LITELLM_CONFIG_FILE")
if not litellm_config_file:
    sys.exit("LITELLM_CONFIG_FILE environment variable is not set.")
with open(litellm_config_file, "r") as f:
    config = yaml.safe_load(f)

servers = config.get("servers", {})
if servers and provider not in servers:
    sys.exit(f"Provider '{provider}' not found in configuration.")

# Select the source: specific provider config if servers exist, otherwise global config (backward compatibility)
source = servers[provider] if servers else config

litellm_api_key = source.get("LITELLM_API_KEY")
litellm_base_url = source.get("LITELLM_BASE_URL")

if not litellm_api_key:
    sys.exit(
        "LiteLLM API key is not configured! Please set LITELLM_API_KEY environment variable."
    )

if not litellm_base_url:
    sys.exit(
        "LiteLLM base URL is not configured! Please set LITELLM_BASE_URL environment variable."
    )

# LLM generation can be slow for large contexts; allow a generous, configurable
# timeout. We disable the SDK's internal retries so the loop below owns backoff
# and logging. Override via LITELLM_REQUEST_TIMEOUT / MAX_RETRIES / MAX_DELAY in
# the config file or environment.
request_timeout = float(os.environ.get("LITELLM_REQUEST_TIMEOUT", "600"))

client = OpenAI(
    api_key=litellm_api_key,
    base_url=litellm_base_url,
    timeout=request_timeout,
    max_retries=0,
)


def read_text_file(file_path):
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return f.read()
    except UnicodeDecodeError:
        try:
            with open(file_path, "r", encoding="latin-1") as f:
                return f.read()
        except Exception:
            sys.exit(f"Could not read file {file_path} as text")


def get_image_mime_type(image_path):
    import mimetypes

    mime_type, _ = mimetypes.guess_type(image_path)
    if mime_type and mime_type.startswith("image/"):
        return mime_type
    if image_path.lower().endswith((".png", ".jpg", ".jpeg", ".gif", ".tiff", ".bmp")):
        ext = image_path.lower().split(".")[-1]
        if ext == "jpg":
            ext = "jpeg"
        return f"image/{ext}"
    return "image/jpeg"


def encode_image_to_base64(image_path):
    import base64

    try:
        with open(image_path, "rb") as image_file:
            base64_image = base64.b64encode(image_file.read()).decode("utf-8")
        mime_type = get_image_mime_type(image_path)
        return f"data:{mime_type};base64,{base64_image}"
    except Exception:
        sys.exit(f"Could not process image file: {image_path}")


valid_model_types = {
    "text": {"text"},
    "image": {"image"},
    "multimodal": {"text", "image"},
}

if model_type not in valid_model_types:
    sys.exit(
        f"Invalid model_type '{model_type}'. Must be one of: {', '.join(valid_model_types)}"
    )

contents = []
for file_path, file_type in context_files:
    if file_type not in valid_model_types[model_type]:
        sys.exit(f"File type '{file_type}' not allowed for model_type '{model_type}'.")
    if file_type == "image":
        contents.append(
            {
                "type": "image_url",
                "image_url": {"url": encode_image_to_base64(file_path)},
            }
        )
    else:
        contents.append(
            {
                "type": "text",
                "text": f"File: {file_path}\nContent:\n{read_text_file(file_path)}",
            }
        )

if question and "text" in valid_model_types[model_type]:
    contents.append({"type": "text", "text": question})

if not contents:
    sys.exit("No input content provided.")

messages = [{"role": "user", "content": contents}]


max_retries = int(config.get("MAX_RETRIES", 3))
max_delay = float(config.get("MAX_DELAY", 900))

# Transient errors that are safe to retry: timeouts and connection failures
# (network/proxy hiccups), rate limiting, and upstream 5xx.
retryable_errors = (
    APITimeoutError,
    APIConnectionError,
    RateLimitError,
    InternalServerError,
)

# An input larger than the model's context window is rejected up front with an
# HTTP 400 -- by the serving backend (vLLM: "Input length (270080) exceeds
# model's maximum context length (128000)"), or by LiteLLM's own pre-call check
# when the proxy sets `enable_pre_call_checks` plus `model_info.max_input_tokens`
# ("Max Input Tokens=..., Got=..."). Either way the rejection is authoritative
# and carries the real numbers, so no local token estimate is needed. The OpenAI
# SDK surfaces it as BadRequestError; without the handler below it escapes the
# retry loop as a raw traceback.
CONTEXT_OVERFLOW_MARKERS = (
    "maximum context length",  # vLLM and OpenAI phrasings
    "context window",  # "context window exceeded" phrasings
    "max input tokens",  # litellm pre-call check
    "too many tokens",
)


def is_context_overflow(exc):
    """True if a 400 is the model's context window being exceeded."""
    message = str(exc).lower()
    return any(marker in message for marker in CONTEXT_OVERFLOW_MARKERS)


# Timeouts on a long generation usually mean the work exceeds the per-request
# budget, so re-sending reproduces the same timeout.  Cap timeout retries
# separately (other transient errors keep the full max_retries budget).
max_timeout_retries = int(os.environ.get("MAX_TIMEOUT_RETRIES") or config.get("MAX_TIMEOUT_RETRIES", 1))
timeout_attempts = 0

# We send no max_tokens: a fixed cap is dangerous for reasoning models (e.g.
# OpenAI o-series, DeepSeek-R1, GLM thinking models), which spend a large,
# unpredictable token budget on chain-of-thought before emitting any answer, so
# a small cap yields finish_reason='length' with empty content.  Rely on the
# model/proxy default instead; truncation is surfaced explicitly below.


def stream_completion():
    """Run one streaming chat completion. Returns (answer_text, finish_reason).

    Streaming keeps the connection alive token-by-token (reasoning models stream
    their chain-of-thought continuously), so request_timeout bounds inter-chunk
    idle rather than total walltime -- the fix for the 600s idle timeout that
    killed long non-streaming jobs.  The answer arrives on `content`;
    `reasoning_content` (chain-of-thought, absent on non-reasoning models) is
    ignored since a Galaxy tool produces a single dataset.
    """
    api_params = {"model": model, "messages": messages, "stream": True}
    # Attribute the request to the Galaxy user so proxies (e.g. LiteLLM) can
    # meter usage and apply per-user budgets/rate limits, via the standard
    # OpenAI `user` field. The id is namespaced by the Galaxy instance URL and
    # hashed: the URL keeps ids unique when several Galaxy instances share one
    # proxy (e.g. usegalaxy.eu), and hashing means no instance-identifying or
    # personal data leaves for the provider. Anonymous users (no id) fall back
    # to a per-instance shared "anonymous" bucket.
    raw_user = galaxy_user_id or "anonymous"
    api_params["user"] = hashlib.sha256(f"{galaxy_url}|{raw_user}".encode()).hexdigest()
    if temperature is not None:
        api_params["temperature"] = temperature

    answer_parts = []
    finish_reason = None
    with client.chat.completions.create(**api_params) as stream:
        for chunk in stream:
            if not chunk.choices:
                continue
            choice = chunk.choices[0]
            delta = choice.delta
            piece = getattr(delta, "content", None)
            if piece:
                answer_parts.append(piece)
            # Discard reasoning_content (chain-of-thought); keep only the answer.
            _ = getattr(delta, "reasoning_content", None)
            if choice.finish_reason is not None:
                finish_reason = choice.finish_reason
    return "".join(answer_parts), finish_reason


for attempt in range(max_retries):
    try:
        answer, finish_reason = stream_completion()

        # An empty answer is never a usable dataset -- fail loudly before
        # writing instead of silently emitting a blank output.md.
        if not answer:
            if finish_reason == "length":
                sys.exit(
                    "The model exhausted its token budget before producing an "
                    "answer (finish_reason='length', empty content). The input "
                    "is likely too large for the model/proxy token limit; split "
                    "it into smaller chunks and re-run."
                )
            sys.exit(
                "The model returned an empty response with no answer content "
                f"(finish_reason={finish_reason!r}). The upstream model/proxy "
                "may have dropped the stream; try again or use a different model."
            )

        # Warn (but still write) when the result may be partial: a token-budget
        # cutoff, upstream content filtering, or a stream that ended without an
        # explicit stop signal.
        if finish_reason == "length":
            print(
                "WARNING: model output was truncated (finish_reason='length'). "
                "The written file is a PARTIAL result. For large inputs, split "
                "them into smaller chunks and re-run.",
                file=sys.stderr,
            )
        elif finish_reason == "content_filter":
            print(
                "WARNING: the model output was filtered by the upstream "
                "model/proxy (finish_reason='content_filter'). The written file "
                "may be partial or altered; please verify and re-run if needed.",
                file=sys.stderr,
            )
        elif finish_reason is None:
            print(
                "WARNING: the stream ended without an explicit stop signal "
                "(finish_reason is None). The written result may be incomplete; "
                "please verify and re-run if needed.",
                file=sys.stderr,
            )

        with open("output.md", "w") as f:
            f.write(answer)
        break
    except BadRequestError as e:
        # The request itself is invalid, so re-sending it reproduces the same
        # rejection -- never retry. Exit with an actionable message instead of
        # letting the SDK exception escape as a traceback.
        if is_context_overflow(e):
            sys.exit(
                "The input is too large for the selected model's context window, "
                "so the request was rejected before any generation started. "
                "Split the input into smaller chunks and re-run, or select a "
                f"model with a larger context window. Upstream error: {e}"
            )
        sys.exit(f"The model/proxy rejected the request: {e}")
    except APITimeoutError as e:
        timeout_attempts += 1
        if attempt == max_retries - 1 or timeout_attempts > max_timeout_retries:
            sys.exit(
                f"Stopped after {timeout_attempts} timeout(s) "
                f"(cap {max_timeout_retries}). Last error: "
                f"{type(e).__name__}: {e}. The request may be too large for the "
                f"timeout budget ({request_timeout}s). Increase "
                f"LITELLM_REQUEST_TIMEOUT or reduce the input size."
            )
        sleep_time = min(2**attempt + random.uniform(0, 1), max_delay)
        print(
            f"{type(e).__name__} encountered ({e}). Timeout attempt "
            f"{timeout_attempts}/{max_timeout_retries}; retrying in "
            f"{sleep_time:.2f}s...",
            file=sys.stderr,
        )
        time.sleep(sleep_time)
    except retryable_errors as e:
        if attempt == max_retries - 1:
            sys.exit(
                f"Max retries ({max_retries}) reached. Last error: "
                f"{type(e).__name__}: {e}"
            )
        sleep_time = min(2**attempt + random.uniform(0, 1), max_delay)
        if isinstance(e, RateLimitError) and hasattr(e, "response") and e.response is not None:
            retry_after = e.response.headers.get("retry-after")
            if retry_after:
                sleep_time = min(float(retry_after), max_delay)
        print(
            f"{type(e).__name__} encountered ({e}). Retrying in "
            f"{sleep_time:.2f} seconds...",
            file=sys.stderr,
        )
        time.sleep(sleep_time)
    except APIError as e:
        # Everything the SDK does not map to a named class handled above: 413
        # (payload too large), 401/403/404/409/422, and APIResponseValidationError
        # (an APIError but NOT an APIStatusError, so a status-based catch misses
        # it). None are worth retrying -- re-sending reproduces them. Must stay
        # LAST: APITimeoutError, APIConnectionError, RateLimitError and
        # InternalServerError are all APIError subclasses caught above.
        sys.exit(f"The request failed and cannot be retried: {type(e).__name__}: {e}")
