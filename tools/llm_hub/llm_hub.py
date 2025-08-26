import json
import os
import random
import sys
import time

import yaml
from openai import InternalServerError, OpenAI

context_files = json.loads(sys.argv[1])
question = sys.argv[2]
model = sys.argv[3]
model_type = sys.argv[4]

litellm_config_file = os.environ.get("LITELLM_CONFIG_FILE")
if not litellm_config_file:
    sys.exit("LITELLM_CONFIG_FILE environment variable is not set.")
with open(litellm_config_file, "r") as f:
    config = yaml.safe_load(f)

litellm_api_key = config.get("LITELLM_API_KEY")
litellm_base_url = config.get("LITELLM_BASE_URL")

if not litellm_api_key:
    sys.exit(
        "LiteLLM API key is not configured! Please set LITELLM_API_KEY environment variable."
    )

if not litellm_base_url:
    sys.exit(
        "LiteLLM base URL is not configured! Please set LITELLM_BASE_URL environment variable."
    )

client = OpenAI(
    api_key=litellm_api_key,
    base_url=litellm_base_url,
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
    if image_path.lower().endswith((".png", ".jpg", ".jpeg", ".gif")):
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


max_retries = config.get("MAX_RETRIES", 3)
max_delay = config.get("MAX_DELAY", 900)
for attempt in range(max_retries):
    try:
        response = client.chat.completions.create(model=model, messages=messages)
        with open("output.md", "w") as f:
            f.write(response.choices[0].message.content or "")
        break
    except InternalServerError as e:
        if attempt == max_retries - 1:
            sys.exit("Max retries reached. Exiting.")
        sleep_time = min(2**attempt + random.uniform(0, 1), max_delay)
        print(
            f"InternalServerError encountered ({e}). Retrying in {sleep_time:.2f} seconds..."
        )
        time.sleep(sleep_time)
