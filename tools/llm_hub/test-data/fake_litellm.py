"""Fake LiteLLM proxy for the tool tests only (Python standard library).

Serves the routes llm_hub uses: /model/info, /utils/token_counter (counts
words) and streaming /chat/completions (answers "OK"). Once listening, it
writes a llm_hub config file pointing to itself, then serves until killed
(or for at most 10 minutes).
"""

import json
import os
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

MODELS = [
    {"model_name": "test-model", "model_info": {"max_input_tokens": 50}},
    {"model_name": "no-limit-model", "model_info": {}},
]


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def send_json(self, data):
        body = json.dumps(data).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.endswith("/model/info"):
            return self.send_json({"data": MODELS})
        self.send_error(404)

    def do_POST(self):
        request = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
        if self.path.endswith("/utils/token_counter"):
            texts = []
            for message in request.get("messages", []):
                content = message.get("content")
                parts = content if isinstance(content, list) else [{"type": "text", "text": content}]
                texts += [p.get("text") or "" for p in parts if p.get("type") == "text"]
            return self.send_json({"total_tokens": len(" ".join(texts).split()), "tokenizer_type": "huggingface_tokenizer"})
        if self.path.endswith("/chat/completions"):
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.end_headers()
            for content, finish in (("OK", None), ("", "stop")):
                chunk = {"id": "fake", "object": "chat.completion.chunk", "created": 0, "model": request.get("model"),
                         "choices": [{"index": 0, "delta": {"content": content}, "finish_reason": finish}]}
                self.wfile.write(f"data: {json.dumps(chunk)}\n\n".encode())
            self.wfile.write(b"data: [DONE]\n\n")
            return
        self.send_error(404)


server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
timer = threading.Timer(600, server.shutdown)
timer.daemon = True
timer.start()
with open(sys.argv[1] + ".tmp", "w") as f:
    f.write(f'LITELLM_API_KEY: test\nLITELLM_BASE_URL: "http://127.0.0.1:{server.server_port}"\n')
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
server.serve_forever()
