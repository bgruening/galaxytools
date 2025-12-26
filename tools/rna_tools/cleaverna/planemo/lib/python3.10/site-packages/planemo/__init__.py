import importlib.metadata

planemo_metadata = importlib.metadata.metadata("planemo")

__version__ = "0.75.33"

PROJECT_NAME = planemo_metadata["Name"]
PROJECT_EMAIL = planemo_metadata["Author-email"].split(" ")[-1]
PROJECT_AUTHOR = PROJECT_USERNAME = "galaxyproject"

PROJECT_URL = "https://github.com/galaxyproject/planemo"
RAW_CONTENT_URL = f"https://raw.github.com/{PROJECT_USERNAME}/{PROJECT_NAME}/master/"
