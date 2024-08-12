import os
import sys

from openai import OpenAI

context_files = sys.argv[1].split(",")
question = sys.argv[2]
model = sys.argv[3]
with open(sys.argv[4], "r") as f:
    openai_api_key = f.read().strip()
if not openai_api_key:
    raise Exception("OpenAI API key is not provided in user preferences!")

client = OpenAI(api_key=openai_api_key)

file_search_sup_ext = [
    "c",
    "cs",
    "cpp",
    "doc",
    "docx",
    "html",
    "java",
    "json",
    "md",
    "pdf",
    "php",
    "pptx",
    "py",
    "rb",
    "tex",
    "txt",
    "css",
    "js",
    "sh",
    "ts",
]

vision_sup_ext = ["jpg", "jpeg", "png", "webp", "gif"]

file_search_file_streams = []
image_urls = []

for path in context_files:
    ext = path.split(".")[-1].lower()
    if ext in vision_sup_ext and model in ["gpt-4o", "gpt-4o-mini", "gpt-4-turbo"]:
        if os.path.getsize(path) > 20 * 1024 * 1024:
            raise Exception(f"File {path} exceeds the 20MB limit and will not be processed.")
        file = client.files.create(file=open(path, "rb"), purpose="vision")
        promt = {"type": "image_file", "image_file": {"file_id": file.id}}
        image_urls.append(promt)

    elif ext in file_search_sup_ext:
        file_search_file_streams.append(open(path, "rb"))
    else:
        raise Exception("Not supported file!")

assistant = client.beta.assistants.create(
    instructions="You are going to get question about the file(s).",
    model=model,
    tools=[{"type": "file_search"}] if file_search_file_streams else None,
)
if file_search_file_streams:
    vector_store = client.beta.vector_stores.create()
    file_batch = client.beta.vector_stores.file_batches.upload_and_poll(
        vector_store_id=vector_store.id, files=file_search_file_streams
    )
    assistant = client.beta.assistants.update(
        assistant_id=assistant.id, tool_resources={"file_search": {"vector_store_ids": [vector_store.id]}}
    )

messages = [
    {
        "role": "user",
        "content": [
            {
                "type": "text",
                "text": question,
            },
            *image_urls,
        ],
    }
]
thread = client.beta.threads.create(messages=messages)
run = client.beta.threads.runs.create_and_poll(thread_id=thread.id, assistant_id=assistant.id)
messages = list(client.beta.threads.messages.list(thread_id=thread.id, run_id=run.id))

message_content = messages[0].content[0].text.value
print("Output has been saved!")
with open("output.txt", "w") as f:
    f.write(message_content)

for image in image_urls:
    client.files.delete(image["image_file"]["file_id"])
if file_search_file_streams:
    client.beta.vector_stores.delete(vector_store.id)
client.beta.threads.delete(thread.id)
client.beta.assistants.delete(assistant.id)
