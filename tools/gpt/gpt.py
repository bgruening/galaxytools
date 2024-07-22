import sys

from openai import OpenAI

context = sys.argv[1]
question = sys.argv[2]
model = sys.argv[3]
openai_api_key = sys.argv[4]
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
