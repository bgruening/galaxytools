import os

import yaml
from bioblend import galaxy

__version__ = "0.10.11"

PROJECT_NAME = "ephemeris"
PROJECT_OWNER = PROJECT_USERAME = "galaxyproject"
PROJECT_URL = "https://github.com/galaxyproject/ephemeris"
PROJECT_AUTHOR = "Galaxy Project and Community"
PROJECT_EMAIL = "jmchilton@gmail.com"
RAW_CONTENT_URL = f"https://raw.github.com/{PROJECT_USERAME}/{PROJECT_NAME}/master/"


def get_or_create_history(history_name: str, gi: galaxy.GalaxyInstance):
    histories = gi.histories.get_histories(name=history_name)
    if histories:
        return histories[0]
    else:
        return gi.histories.create_history(name=history_name)


def check_url(url, log=None):
    if not url.startswith("http"):
        if log:
            log.warning("URL should start with http:// or https://. https:// chosen by default.")
        url = "https://" + url
    return url


def get_galaxy_connection(args, file=None, log=None, login_required=True):
    """
    Return a Galaxy connection, given a user or an API key.
    If not given gets the arguments from the file.
    If either is missing raise ValueError.
    """
    if file:
        file_content = load_yaml_file(file)
    else:
        file_content = dict()

    url = args.galaxy or file_content.get("galaxy_instance")
    galaxy_url = check_url(url, log)
    api_key = args.api_key or file_content.get("api_key") or os.environ.get("EPHEMERIS_API_KEY")

    if args.user and args.password:
        return galaxy.GalaxyInstance(url=galaxy_url, email=args.user, password=args.password)
    elif api_key:
        return galaxy.GalaxyInstance(url=galaxy_url, key=api_key)
    elif not login_required:
        return galaxy.GalaxyInstance(url=galaxy_url)
    else:
        raise ValueError("Missing api key or user & password combination, in order to make a galaxy connection.")


def load_yaml_file(filename):
    """
    Load YAML from the `tool_list_file` and return a dict with the content.
    """
    with open(filename) as f:
        dictionary = yaml.safe_load(f)
    return dictionary


def dump_to_yaml_file(content, file_name):
    """
    Dump YAML-compatible `content` to `file_name`.
    """
    with open(file_name, "w") as f:
        yaml.dump(content, f, default_flow_style=False)
