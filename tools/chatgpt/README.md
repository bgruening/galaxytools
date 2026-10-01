# chatGPT Galaxy tool

Sends your prompt, and optional files, to a large language model (LLM) and saves the answer as a Markdown file.
You can use OpenAI, or any server with an OpenAI-compatible API, for example Open WebUI, vLLM, Ollama or LiteLLM.

## Credentials

Add them once with the **Provide credentials** button at the top of the tool form.

- **OpenAI**: your OpenAI API key ([get one here](https://platform.openai.com/account/api-keys); the account needs [credits](https://platform.openai.com/settings/organization/billing)).
- **Custom server**: the server URL and, if the server needs one, an API key. The URL is the server address plus its API path, for example `https://my-server.org/v1`. For Open WebUI it ends with `/api`. The job runs on a Galaxy server, so `localhost` means that server, not your computer.

## Inputs

- **Server**: OpenAI or your custom server.
- **Model**: for OpenAI, choose from the current models or type any model name. The list is loaded live from the public [OpenRouter model list](https://openrouter.ai/api/v1/models), because OpenAI's own model list needs the user's key, which Galaxy cannot use while building the form. The `*-pro` and `*-codex` models are left out: they only work with OpenAI's Responses API. For a custom server, type the model name as the server knows it.
- **Context** (optional): files for the model to read. Text files (TXT, CSV, JSON, HTML) work best. PDF and Word files are not converted, so turn them into text first, for example with the Markitdown tool. Images (JPG, PNG, GIF, max 20 MB each) need a model that can read images: some models ignore images without an error.
- **Prompt**: your question or task. Be specific.

## Advanced options

- **Temperature**: lower (0 to 0.3) gives focused, repeatable answers. Higher (0.7 or more) gives more creative answers.
- **Top P**: another way to control randomness. Change Temperature or Top P, not both.
- **Max tokens**: the longest answer allowed. Reasoning models also count their hidden thinking, so a low value can give an empty answer.
- **System message**: tells the model how to behave, for example "You are a helpful biology assistant".

Many reasoning models do not accept Temperature or Top P. The tool then leaves them out and writes a note in the job log.

## If something goes wrong

The job log says what went wrong in one short sentence, for example a wrong API key, an unknown model name, or files that are too long for the model.

## Privacy

Your prompt, system message and files leave Galaxy. They are sent to OpenAI or to your custom server, and that service's data policy applies. Files are sent inside the request (text inline, images as base64), not uploaded to OpenAI's file storage.
