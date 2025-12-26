"""Module contains code for the Topic class, dealing with the creation of a training topic."""

import collections
import os

from galaxy.util import listify

from planemo import templates
from .utils import (
    load_yaml,
    Requirement,
    save_to_yaml,
)

INDEX_FILE_TEMPLATE = """---
layout: topic
topic_name: {{ topic }}
---
"""

README_FILE_TEMPLATE = """
{{ topic }}
==========

Please refer to the [CONTRIBUTING.md](../../CONTRIBUTING.md) before adding or updating any material
"""


DOCKER_FILE_TEMPLATE = """
# Galaxy - {{ topic_title }}
#
# to build the docker image, go to root of training repo and
#    docker build -t {{ topic_name }} -f topics/{{ topic_name }}/docker/Dockerfile .
#
# to run image:
#    docker run -p "8080:80" -t {{ topic_name }}
#    use -d to automatically dowload the datalibraries in the container

FROM quay.io/bgruening/galaxy:20.05

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: {{ topic_title }}"

# copy the tutorials directory for your topic
ADD topics/{{ topic_name }}/tutorials/ /tutorials/

# install everything for tutorials
ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
ADD bin/mergeyaml.py /mergeyaml.py
ADD bin/data_libarary_download.sh /data_libarary_download.sh
RUN /setup-tutorials.sh

ENTRYPOINT ["/data_libarary_download.sh"]
"""


class Topic:
    """Class to describe a training topic."""

    def __init__(self, name="new_topic", target="use", title="The new topic", summary="Summary", parent_dir="topics"):
        """Init a topic instance."""
        self.name = name
        self.type = target
        self.title = title
        self.summary = summary
        self.docker_image = ""
        self.editorial_board = []
        self.parent_dir = parent_dir
        self.set_default_requirement()
        self.set_paths()

    def init_from_kwds(self, kwds):
        """Init a topic instance from a kwds dictionary."""
        self.name = kwds["topic_name"]
        self.type = kwds["topic_target"]
        self.title = kwds["topic_title"]
        self.summary = kwds["topic_summary"]
        self.set_default_requirement()
        self.set_paths()

    def init_from_metadata(self):
        """Init a topic instance from the metadata file."""
        metadata = load_yaml(self.metadata_fp)
        self.name = metadata["name"]
        self.type = metadata["type"]
        self.title = metadata["title"]
        self.summary = metadata["summary"]
        self.requirements = []
        if "requirements" in metadata:
            for r in listify(metadata["requirements"]):
                req = Requirement()
                req.init_from_dict(r)
                self.requirements.append(req)
        if "docker_image" in metadata:
            self.docker_image = metadata["docker_image"]
        self.editorial_board = metadata["editorial_board"]
        self.set_paths()

    # GETTERS
    def get_requirements(self):
        """Get the requirements as a list of ordered dictionaries."""
        reqs = []
        for req in self.requirements:
            reqs.append(req.export_to_ordered_dict())
        return reqs

    def export_metadata_to_ordered_dict(self):
        """Export the topic metadata into an ordered dictionary."""
        metadata = collections.OrderedDict()
        metadata["name"] = self.name
        metadata["type"] = self.type
        metadata["title"] = self.title
        metadata["summary"] = self.summary
        metadata["requirements"] = self.get_requirements()
        metadata["docker_image"] = self.docker_image
        metadata["editorial_board"] = self.editorial_board
        return metadata

    # SETTERS
    def set_default_requirement(self):
        """Set default requirement: Galaxy introduction."""
        self.requirements = []
        if self.type == "use":
            self.requirements.append(Requirement())

    def set_paths(self):
        """Set the paths to folder and files."""
        self.dir = os.path.join(self.parent_dir, self.name)
        self.img_folder = os.path.join(self.dir, "images")
        self.tuto_folder = os.path.join(self.dir, "tutorials")
        self.index_fp = os.path.join(self.dir, "index.md")
        self.readme_fp = os.path.join(self.dir, "README.md")
        self.metadata_fp = os.path.join(self.dir, "metadata.yaml")
        self.docker_folder = os.path.join(self.dir, "docker")
        self.dockerfile_fp = os.path.join(self.docker_folder, "Dockerfile")

    # TESTS
    def exists(self):
        """Test if the topic exists."""
        return os.path.isdir(self.dir)

    # OTHER METHODS
    def create_topic_structure(self):
        """Create the skeleton of a new topic.

        1. create the folder and its structure
        2. update the index.md to match your topic's name
        3. fill the metadata
        4. add a symbolic link to the metadata.yaml from the metadata folder
        """
        # create the folder and its structure
        os.makedirs(self.dir)
        self.img_folder = os.path.join(self.dir, "images")
        os.makedirs(self.img_folder)
        self.tuto_folder = os.path.join(self.dir, "tutorials")
        os.makedirs(self.tuto_folder)

        # create the index.md and add the topic name
        self.index_fp = os.path.join(self.dir, "index.md")
        with open(self.index_fp, "w") as index_f:
            index_f.write(templates.render(INDEX_FILE_TEMPLATE, **{"topic": self.name}))

        # create the README file
        self.readme_fp = os.path.join(self.dir, "README.md")
        with open(self.readme_fp, "w") as readme_f:
            readme_f.write(templates.render(README_FILE_TEMPLATE, **{"topic": self.title}))

        # create the metadata file
        self.metadata_fp = os.path.join(self.dir, "metadata.yaml")
        save_to_yaml(self.export_metadata_to_ordered_dict(), self.metadata_fp)

        # create Dockerfile
        self.docker_folder = os.path.join(self.dir, "docker")
        os.makedirs(self.docker_folder)
        self.dockerfile_fp = os.path.join(self.docker_folder, "Dockerfile")
        with open(self.dockerfile_fp, "w") as dockerfile:
            dockerfile.write(
                templates.render(DOCKER_FILE_TEMPLATE, **{"topic_name": self.name, "topic_title": self.title})
            )

        # add a symbolic link to the metadata.yaml
        metadata_dir = "metadata"
        if not os.path.isdir(metadata_dir):
            os.makedirs(metadata_dir)
        os.chdir(metadata_dir)
        os.symlink(os.path.join("..", self.metadata_fp), "%s.yaml" % self.name)
        os.chdir("..")
