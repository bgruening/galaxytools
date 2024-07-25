import sys

from openai import OpenAI

context = sys.argv[1]
question = sys.argv[2]
model = sys.argv[3]
openai_api_key_path = sys.argv[4]
with open(openai_api_key_path, "r") as f:
    openai_api_key = f.read().strip()
client = OpenAI(api_key=openai_api_key)
response = client.chat.completions.create(
    model=model,
    messages=[
        {
            "role": "system",
            "content": f"You are going to get question about this text: {context}",
        },
        {"role": "user", "content": f"{question}"},
    ],
)

print(response.choices[0].message.content)
