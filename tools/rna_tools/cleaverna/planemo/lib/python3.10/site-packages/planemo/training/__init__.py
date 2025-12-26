"""Module contains code for gtdk: Galaxy training development kit."""

from planemo.io import info
from .topic import Topic
from .tutorial import Tutorial


class Training:
    """Class to describe a training."""

    def __init__(self, kwds):
        """Init an instance of Training."""
        self.kwds = kwds
        self.topics_dir = "topics"
        self.topic = Topic(parent_dir=self.topics_dir, name=kwds["topic_name"])
        self.galaxy_url = kwds["galaxy_url"] if "galaxy_url" in kwds else ""
        self.galaxy_api_key = kwds["galaxy_api_key"] if "galaxy_api_key" in kwds else ""
        self.tuto = None

    def init_training(self, ctx):
        """Create/update a topic/tutorial."""
        if not self.topic.exists():
            info("The topic %s does not exist. It will be created" % self.topic.name)
            self.topic.init_from_kwds(self.kwds)
            self.topic.create_topic_structure()

        if not self.kwds["tutorial_name"]:
            if self.kwds["slides"]:
                raise Exception("A tutorial name is needed to create the skeleton of a tutorial slide deck")
            if self.kwds["workflow"] or self.kwds["workflow_id"]:
                raise Exception("A tutorial name is needed to create the skeleton of the tutorial from a workflow")
            if self.kwds["zenodo_link"]:
                raise Exception("A tutorial name is needed to add Zenodo information")
        else:
            self.tuto = Tutorial(training=self, topic=self.topic)
            self.tuto.init_from_kwds(self.kwds)
            if not self.tuto.exists():
                info(f"The tutorial {self.tuto.name} in topic {self.topic.name} does not exist. It will be created.")
                self.tuto.create_tutorial(ctx)
        info(
            "WARNING: Change the contributors listed in the metadata of the new training before serving the website to fit the one listed in the CONTRIBUTORS.yaml file"
        )

    def check_topic_init_tuto(self):
        """Check that the topic and tutorial are already there and retrieve them."""
        # check topic
        if not self.topic.exists():
            raise Exception("The topic %s does not exists. It should be created" % self.topic.name)
        self.topic.init_from_metadata()
        # initiate the tutorial
        self.tuto = Tutorial(training=self, topic=self.topic)
        self.tuto.init_from_existing_tutorial(self.kwds["tutorial_name"])
        if "workflow" in self.kwds:
            self.tuto.init_wf_fp = self.kwds["workflow"]
        if "workflow_id" in self.kwds:
            self.tuto.init_wf_id = self.kwds["workflow_id"]

    def fill_data_library(self, ctx):
        """Fill a data library for a tutorial."""
        self.check_topic_init_tuto()
        # get the zenodo link
        z_link = ""
        update_metadata = False
        if self.tuto.zenodo_link != "":
            if self.kwds["zenodo_link"]:
                info("The data library and the metadata will be updated with the new Zenodo link")
                z_link = self.kwds["zenodo_link"]
                self.tuto.zenodo_link = z_link
                update_metadata = True
            else:
                info("The data library will be extracted using the Zenodo link in the metadata of the tutorial")
                z_link = self.tuto.zenodo_link
        elif self.kwds["zenodo_link"]:
            info("The data library will be created and the metadata will be filled with the new Zenodo link")
            z_link = self.kwds["zenodo_link"]
            self.tuto.zenodo_link = z_link
            update_metadata = True

        if z_link == "" or z_link is None:
            raise Exception(
                "A Zenodo link should be provided either in the metadata file or as argument of the command"
            )

        # extract the data library from Zenodo
        self.tuto.prepare_data_library_from_zenodo()

        # update the metadata
        if update_metadata:
            self.tuto.write_hands_on_tutorial(add_z_file_links=False)

    def generate_tuto_from_wf(self, ctx):
        """Generate the skeleton of a tutorial from a workflow."""
        self.check_topic_init_tuto()
        if self.tuto.has_workflow():
            info("Create tutorial skeleton from workflow")
            self.tuto.create_hands_on_tutorial(ctx)
            self.tuto.export_workflow_file()
        else:
            raise Exception(
                "A path to a local workflow or the id of a workflow on a running Galaxy instance should be provided"
            )
