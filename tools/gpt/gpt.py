import sys

from openai import OpenAI

context = sys.argv[1]
question = sys.argv[2]
client = OpenAI(api_key="")

response = client.chat.completions.create(
    model="gpt-4o",
    messages=[
        {
            "role": "system",
            "content": f"You are going to get question about this text: {context}",
        },
        {"role": "user", "content": f"{question}"},
    ],
)

print(response.choices[0].message.content)
