# LangChain Text Splitters Galaxy Tool Dev Utilities

## Populate tiktoken cache directory

Download the OpenAI encoding files using the script `get_encodings.py` in this folder.
This will populate the cache directory `/tmp/tiktoken-cache`.

## Testing tiktoken without internet access in Galaxy

! Prerequisite: The OpenAI encoding files must be downloaded and the cache directory populated as described above.

To test whether a Galaxy tool can use a locally cached tiktoken encoding without internet access, 
save the following configuration as job_conf_offline.xml:

```xml
<?xml version="1.0"?>
<job_conf>
    <plugins>
        <plugin 
            id="local"
            type="runner"
            load="galaxy.jobs.runners.local:LocalJobRunner"
            workers="4"
        />
    </plugins>

    <destinations default="local">
        <destination id="local" runner="local">
            <env id="TIKTOKEN_CACHE_DIR">/tmp/tiktoken-cache</env>
            <!-- Point to an unreachable local port -->
            <env id="HTTP_PROXY">http://127.0.0.1:9</env>
            <env id="HTTPS_PROXY">http://127.0.0.1:9</env>
        </destination>
    </destinations>
</job_conf>
```
Start Planemo using the job configuration above:
`planemo serve --job_config_file dev_utils/job_conf_offline.xml`

You can test as follows that the tool works without internet access by setting
**Text splitting boundaries** to:

**(A) Between tokens**.

**(B) Character**, then set **Target chunk size should be determined by the count of** to **Tokens**.

To double check that the tool actually fails when the cache is empty, remove or rename the cache directory `/tmp/tiktoken-cache` and run it as described above. It should run into a 
`ConnectionRefusedError` when trying to download the encoding files from the OpenAI servers.