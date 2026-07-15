import json
import os
import random
import sys
import time

import yaml
from openai import (
    APIConnectionError,
    APITimeoutError,
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
