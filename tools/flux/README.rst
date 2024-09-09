FLUX
====

Set up user credentials on Galaxy
---------------------------------

To enable users to set their credentials for this tool, make sure the
file ``config/user_preferences_extra.yml`` has the following section:

::

        preferences:
            flux:
                description: Your FLUX settings
                inputs:
                    - name: huggingface_hub_token
                    label: API Key
                    type: password
                    required: False
