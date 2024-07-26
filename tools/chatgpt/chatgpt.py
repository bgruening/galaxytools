import sys

from openai import OpenAI

context_file = sys.argv[1]
question = sys.argv[2]
model = sys.argv[3]

client = OpenAI()
assistant = client.beta.assistants.create(
    instructions="You are going to get question about the file.", model=model, tools=[{"type": "file_search"}]
)
vector_store = client.beta.vector_stores.create()
file_streams = [open(path, "rb") for path in context_file.split(",")]
file_batch = client.beta.vector_stores.file_batches.upload_and_poll(vector_store_id=vector_store.id, files=file_streams)
assistant = client.beta.assistants.update(
    assistant_id=assistant.id, tool_resources={"file_search": {"vector_store_ids": [vector_store.id]}}
)
thread = client.beta.threads.create(messages=[{"role": "user", "content": question}])
run = client.beta.threads.runs.create_and_poll(thread_id=thread.id, assistant_id=assistant.id)
messages = list(client.beta.threads.messages.list(thread_id=thread.id, run_id=run.id))
message_content = messages[0].content[0].text

response = message_content.value
print("Output has been saved!")
with open("output.txt", "w") as f:
    f.write(response)
