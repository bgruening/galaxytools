import json
import os
import sys

import yaml
from openai import OpenAI

context_files = json.loads(sys.argv[1])
question = sys.argv[2]
model = sys.argv[3]
model_type = sys.argv[4]

litellm_config_file = os.environ.get("LITELLM_CONFIG_FILE")
if not litellm_config_file:
    print("LITELLM_CONFIG_FILE environment variable is not set.")
    sys.exit(1)
with open(litellm_config_file, "r") as f:
    config = yaml.safe_load(f)

litellm_api_key = config.get("LITELLM_API_KEY")
litellm_base_url = config.get("LITELLM_BASE_URL")

if not litellm_api_key:
    print(
        "LiteLLM API key is not configured! Please set LITELLM_API_KEY environment variable."
    )
    sys.exit(1)

if not litellm_base_url:
    print(
        "LiteLLM base URL is not configured! Please set LITELLM_BASE_URL environment variable."
    )
    sys.exit(1)

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
            print(f"Could not read file {file_path} as text")
            sys.exit(1)


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
        print(f"Could not process image file: {image_path}")
        sys.exit(1)


content_text = question
messages = []

valid_model_types = ["text", "multimodal"]
if model_type not in valid_model_types:
    print(
        f"Invalid model_type '{model_type}'. Must be one of: {', '.join(valid_model_types)}"
    )
    sys.exit(1)

if context_files:
    context_text_parts = []
    image_contents = []

    for file_path, file_type in context_files:
        if file_type == "image":
            if model_type == "multimodal":
                base64_image_url = encode_image_to_base64(file_path)
                image_contents.append(
                    {"type": "image_url", "image_url": {"url": base64_image_url}}
                )
            else:
                print(
                    f"Image file '{file_path}' provided, but model_type is not 'multimodal'."
                )
                sys.exit(1)
        else:
            text_content = read_text_file(file_path)
            context_text_parts.append(
                f"File: {file_path}\nContent:\n{text_content}\n---\n"
            )

    if context_text_parts:
        context_text = "Context files:\n\n" + "\n".join(context_text_parts)
        content_text = f"{context_text}\n\nUser Question: {question}"

    if model_type == "multimodal" and image_contents:
        content = [{"type": "text", "text": content_text}, *image_contents]
        messages = [{"role": "user", "content": content}]
    else:
        messages = [{"role": "user", "content": content_text}]
else:
    messages = [{"role": "user", "content": content_text}]

response = client.chat.completions.create(model=model, messages=messages)

with open("output.md", "w") as f:
    f.write(response.choices[0].message.content or "")
