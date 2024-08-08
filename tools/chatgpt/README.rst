ChatGPT
======

Set up user credentials on Galaxy
---------------------------------

To enable users to set their credentials for this tool, make sure the
file ``config/user_preferences_extra.yml`` has the following section:

::

        preferences:
            chatgpt:
                description: Your ChatGPT API settings
                inputs:
                - name: api_key
                    label: API Key
                    type: password
                    required: False
