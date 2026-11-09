"""Fake LiteLLM proxy for the llm_hub speech to text and decision tests only.

Standard library only, so it runs in the tool's own container. It answers
POST /audio/transcriptions like the real proxy: it checks the API key, the
model and the uploaded file, and returns a transcript that names the file,
its size and the language field, so a test can check what reached the server.
With response_format verbose_json it adds two timestamped segments, like
Whisper; test-transcriber-no-timestamps rejects that like Voxtral on vLLM.

POST /systemone answers each question in the System One shape: a yes/no
question is "yes" only if the state contains the word "degraded", a choice
question picks its first option and a score question returns 0.25.

Once listening, it writes a llm_hub config file pointing to itself. It stops
when killed, or after 10 minutes.
"""

import json
import os
import sys
import threading
from email.parser import BytesParser
from email.policy import HTTP
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

API_KEY = "test-key"
MODELS = {"test-transcriber", "test-transcriber-no-timestamps"}
DECISION_MODELS = {"test-decider"}
# like a chat model on a proxy without the System One route
NO_DECISION_MODELS = {"test-transcriber"}


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

    def do_POST(self):
        if self.headers.get("Authorization") != f"Bearer {API_KEY}":
            return self.send_json({"error": {"message": "Authentication Error"}}, 401)
        if self.path.endswith("/systemone"):
            return self.decide()
        if not self.path.endswith("/audio/transcriptions"):
            return self.send_json({"error": {"message": "Not Found"}}, 404)
        body = self.rfile.read(int(self.headers.get("Content-Length", 0)))
        # multipart/form-data, as the OpenAI client sends it
        form = BytesParser(policy=HTTP).parsebytes(
            f"Content-Type: {self.headers.get('Content-Type')}\r\n\r\n".encode() + body
        )
        fields, upload = {}, None
        for part in form.iter_parts():
            if part.get_filename():
                upload = (part.get_filename(), part.get_payload(decode=True))
            else:
                fields[part.get_param("name", header="content-disposition")] = part.get_content().strip()
        if fields.get("model") not in MODELS:
            return self.send_json({"error": {"message": f"Invalid model name passed in model={fields.get('model')}."}}, 400)
        if not upload or not upload[1]:
            return self.send_json({"error": {"message": "No audio file was uploaded."}}, 400)
        if not fields.get("user"):
            return self.send_json({"error": {"message": "No user field was sent."}}, 400)
        # like Whisper and Voxtral, the text starts with a space
        text = f" Transcript of {upload[0]} ({len(upload[1])} bytes, language {fields.get('language') or 'detected'})."
        if fields.get("response_format") != "verbose_json":
            return self.send_json({"text": text})
        if fields["model"] == "test-transcriber-no-timestamps":
            return self.send_json({"error": {"message": "litellm.BadRequestError: Hosted_vllmException - Currently do not support verbose_json for test-transcriber-no-timestamps"}}, 400)
        segments = [{"id": 0, "start": 0.0, "end": 0.5, "text": text}, {"id": 1, "start": 0.5, "end": 1.0, "text": " Second segment."}]
        self.send_json({"text": text, "language": fields.get("language") or "en", "duration": 1.0, "segments": segments})

    def decide(self):
        body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))))
        if body.get("model") in NO_DECISION_MODELS:
            return self.send_json({"error": {"message": "litellm.NotFoundError: NotFoundError: Hosted_vllmException - Not Found"}}, 404)
        if body.get("model") not in DECISION_MODELS:
            return self.send_json({"error": {"message": f"Invalid model name passed in model={body.get('model')}."}}, 400)
        if not body.get("user"):
            return self.send_json({"error": {"message": "No user field was sent."}}, 400)
        state = json.dumps(body.get("state"))
        answers = {}
        for name, q in body.get("questions", {}).items():
            if q["type"] == "noul":
                yes = 1.0 if "degraded" in state else 0.0
                answers[name] = {"type": "noul", "noul": yes, "probabilities": {"yes": yes, "no": 1 - yes}, "confidence": max(yes, 1 - yes)}
            elif q["type"] == "choice":
                labels = list(q["criteria"])
                probabilities = {label: (1.0 if i == 0 else 0.0) for i, label in enumerate(labels)}
                answers[name] = {"type": "choice", "choice": labels[0], "confidence": 1.0, "probabilities": probabilities}
            elif q["type"] == "score":
                legend = {"0": q["criteria"][0], "1": q["criteria"][-1]}
                answers[name] = {"type": "score", "score": 0.25, "confidence": 0.75, "legend": legend, "probabilities": {"0": 0.75, "1": 0.25}}
            else:
                return self.send_json({"error": {"message": f"Unknown question type {q['type']!r}."}}, 400)
        self.send_json({"model": body["model"], "answers": answers, "usage": {"input_tokens": len(state), "output_tokens": 0}, "object": "structured_decision"})


server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
timer = threading.Timer(600, server.shutdown)
timer.daemon = True
timer.start()
with open(sys.argv[1] + ".tmp", "w") as f:
    f.write(f'LITELLM_API_KEY: {API_KEY}\nLITELLM_BASE_URL: "http://127.0.0.1:{server.server_port}"\n')
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
server.serve_forever()
