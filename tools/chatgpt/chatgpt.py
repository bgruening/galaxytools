import ast
import os
import sys

from openai import AuthenticationError, OpenAI

context_files = ast.literal_eval(sys.argv[1])
question = sys.argv[2]
model = sys.argv[3]
with open(sys.argv[4], "r") as f:
    openai_api_key = f.read().strip()
if not openai_api_key:
    print("OpenAI API key is not provided in user preferences!")
    sys.exit(1)

client = OpenAI(api_key=openai_api_key)

file_search_file_streams = []
image_files = []

for path, type in context_files:
    if type == "image":
        if os.path.getsize(path) > 20 * 1024 * 1024:
            print(f"File {path} exceeds the 20MB limit and will not be processed.")
            sys.exit(1)
        file = client.files.create(file=open(path, "rb"), purpose="vision")
        promt = {"type": "image_file", "image_file": {"file_id": file.id}}
        image_files.append(promt)
    else:
        file_search_file_streams.append(open(path, "rb"))

try:
    assistant = client.beta.assistants.create(
        instructions=(
            "You will receive questions about files from file searches "
            "and image files. For file search queries, identify and "
            "retrieve the relevant files based on the question. "
            "For image file queries, analyze the image content and "
            "provide relevant information or insights based on the image data."
        ),
        model=model,
        tools=[{"type": "file_search"}] if file_search_file_streams else [],
    )
except AuthenticationError as e:
    print(f"Authentication error: {e.message}")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {str(e)}")
    sys.exit(1)

if file_search_file_streams:
    vector_store = client.beta.vector_stores.create()
    file_batch = client.beta.vector_stores.file_batches.upload_and_poll(
        vector_store_id=vector_store.id, files=file_search_file_streams
    )
    assistant = client.beta.assistants.update(
        assistant_id=assistant.id,
        tool_resources={"file_search": {"vector_store_ids": [vector_store.id]}},
    )

messages = [
    {
        "role": "user",
        "content": [
            {
                "type": "text",
                "text": question,
            },
            *image_files,
        ],
    }
]
thread = client.beta.threads.create(messages=messages)
run = client.beta.threads.runs.create_and_poll(
    thread_id=thread.id, assistant_id=assistant.id
)
assistant_messages = list(
    client.beta.threads.messages.list(thread_id=thread.id, run_id=run.id)
)
if not assistant_messages:
    print(
        "No output was generated!\nPlease ensure that your OpenAI account has sufficient credits.\n"
        "You can check your balance here: https://platform.openai.com/settings/organization/billing"
    )
    sys.exit(1)
message_content = assistant_messages[0].content[0].text.value
print("Output has been saved!")
with open("output.txt", "w") as f:
    f.write(message_content)

for image in image_files:
    client.files.delete(image["image_file"]["file_id"])
if file_search_file_streams:
    client.beta.vector_stores.delete(vector_store.id)
client.beta.threads.delete(thread.id)
client.beta.assistants.delete(assistant.id)
