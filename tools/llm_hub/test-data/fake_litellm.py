"""Fake LiteLLM proxy in front of vLLM, for the llm_hub tests only.

Standard library only, so it runs in the tool's own container. It behaves
like the real proxy for the routes llm_hub uses, with the same formats:

- GET  /model/info: models with max_input_tokens (when the admin set one).
- POST /utils/token_counter: counts the message content with the model's
  tokenizer (here a simple word/punctuation tokenizer; images count 85).
- POST /chat/completions: checks the API key and the model, rejects images
  for text-only models and prompts over the context length with vLLM's
  HTTP 400 errors, and streams the answer. Like vLLM, the prompt also counts
  the chat template, so it is a few tokens longer than the token counter says.

Once listening, it writes a llm_hub config file pointing to itself. It stops
when killed, or after 10 minutes.
"""

import json
import os
import re
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

API_KEY = "test-key"
TEMPLATE_TOKENS = 5  # chat template tokens vLLM adds around the messages
IMAGE_TOKENS = 85
MODELS = {
    # name: context length of the backend, published limit, accepts images
    "test-model": {"context": 50, "max_input_tokens": 50, "vision": False},
    "no-limit-model": {"context": 50, "max_input_tokens": None, "vision": False},
    "text-only-model": {"context": 500, "max_input_tokens": 500, "vision": False},
    "flaky-model": {"context": 500, "max_input_tokens": None, "vision": False},
}
failures_left = {"flaky-model": 1}  # answer with HTTP 503 this many times first


def count_tokens(messages):
    """Content tokens and number of images in the messages."""
    tokens, images = 0, 0
    for message in messages:
        content = message.get("content")
        parts = content if isinstance(content, list) else [{"type": "text", "text": content or ""}]
        for part in parts:
            if part.get("type") == "image_url":
                images += 1
            else:
                tokens += len(re.findall(r"\w+|[^\w\s]", part.get("text") or ""))
    return tokens + images * IMAGE_TOKENS, images


def litellm_error(message, model):
    return f"litellm.BadRequestError: Hosted_vllmException - {json.dumps({'error': {'message': message, 'type': 'BadRequestError', 'param': None, 'code': 400}})}. Received Model Group={model}"


class Handler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

    def send_json(self, data, status=200):
        body = json.dumps(data).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_error_json(self, message, status):
        self.send_json({"error": {"message": message, "type": None, "param": None, "code": str(status)}}, status)

    def authorized(self):
        if self.headers.get("Authorization") == f"Bearer {API_KEY}":
            return True
        self.send_error_json("Authentication Error, Invalid proxy server token passed.", 401)
        return False

    def do_GET(self):
        if not self.authorized():
            return
        if self.path.endswith("/model/info"):
            data = [{"model_name": name, "litellm_params": {"model": f"hosted_vllm/{name}"},
                     "model_info": {"id": name, "max_input_tokens": spec["max_input_tokens"]}}
                    for name, spec in MODELS.items()]
            return self.send_json({"data": data})
        self.send_error_json("Not Found", 404)

    def do_POST(self):
        if not self.authorized():
            return
        request = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
        model = request.get("model")
        if model not in MODELS:
            return self.send_error_json(f"Invalid model name passed in model={model}. Call `/v1/models` to view available models for your key.", 400)
        spec = MODELS[model]
        tokens, images = count_tokens(request.get("messages", []))
        if self.path.endswith("/utils/token_counter"):
            return self.send_json({"total_tokens": tokens, "request_model": model, "model_used": model,
                                   "tokenizer_type": "huggingface_tokenizer"})
        if not self.path.endswith("/chat/completions"):
            return self.send_error_json("Not Found", 404)
        if images and not spec["vision"]:
            return self.send_error_json(litellm_error(f"{model} is not a multimodal model", model), 400)
        prompt = tokens + TEMPLATE_TOKENS
        if prompt > spec["context"]:
            return self.send_error_json(litellm_error(
                f"This model's maximum context length is {spec['context']} tokens. However, you requested 0 output tokens "
                f"and your prompt contains at least {spec['context'] + 1} input tokens, for a total of at least "
                f"{spec['context'] + 1} tokens. Please reduce the length of the input prompt or the number of requested "
                f"output tokens. (parameter=input_tokens, value={spec['context'] + 1})", model), 400)
        if failures_left.get(model):
            failures_left[model] -= 1
            return self.send_error_json("litellm.InternalServerError: Hosted_vllmException - Service Unavailable", 503)
        answer = {"id": "chatcmpl-fake", "created": 0, "model": model}
        if not request.get("stream"):
            return self.send_json({**answer, "object": "chat.completion",
                                   "choices": [{"index": 0, "message": {"role": "assistant", "content": "OK"}, "finish_reason": "stop"}],
                                   "usage": {"prompt_tokens": prompt, "completion_tokens": 1, "total_tokens": prompt + 1}})
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.end_headers()
        for delta, finish in (({"role": "assistant", "content": "OK"}, None), ({}, "stop")):
            chunk = {**answer, "object": "chat.completion.chunk", "choices": [{"index": 0, "delta": delta, "finish_reason": finish}]}
            self.wfile.write(f"data: {json.dumps(chunk)}\n\n".encode())
        self.wfile.write(b"data: [DONE]\n\n")


server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
timer = threading.Timer(600, server.shutdown)
timer.daemon = True
timer.start()
with open(sys.argv[1] + ".tmp", "w") as f:
    f.write(f'LITELLM_API_KEY: {API_KEY}\nLITELLM_BASE_URL: "http://127.0.0.1:{server.server_port}"\n')
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
server.serve_forever()
