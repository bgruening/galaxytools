"""Module contains code for the Tutorial class, dealing with the creation of a training tutorial."""

import collections
import json
import os
import re
import shutil

import oyaml as yaml
import requests
from bioblend.galaxy import GalaxyInstance

from planemo import templates
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.io import (
    error,
    info,
)
from planemo.runnable import for_path
from .tool_input import (
    get_empty_input,
    get_empty_param,
    ToolInput,
)
from .utils import (
    load_yaml,
    save_to_yaml,
)

TUTO_HAND_ON_TEMPLATE = """---
layout: tutorial_hands_on

{{ metadata }}
---

{{ body }}
"""

TUTO_SLIDES_TEMPLATE = """---
layout: tutorial_slides
logo: "GTN"

{{ metadata }}
---

### How to fill the slide decks?

Please follow our
[tutorial to learn how to fill the slides]({{ '{{' }} site.baseurl {{ '}}' }}/topics/contributing/tutorials/create-new-tutorial-slides/slides.html)
"""


HANDS_ON_TOOL_BOX_TEMPLATE = """
## Sub-step with **{{tool_name}}**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {{ '{%' }} tool [{{tool_name}}]({{tool_id}}) {{ '%}' }} with the following parameters:{{inputlist}}{{paramlist}}
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

"""

TUTO_BIBLIOGRAPHY_TEMPLATE = """
# This is the bibliography file for your tutorial.
#
# To add bibliography (bibtex) entries here, follow these steps:
#  1) Find the DOI for the article you want to cite
#  2) Go to https://doi2bib.org and fill in the DOI
#  3) Copy the resulting bibtex entry into this file
#
# To cite the example below, in your tutorial.md file
# use {{ '{%' }} cite Batut2018 {{ '%}' }}
#
# If you want to cite an online resourse (website etc)
# you can use the 'online' format (see below)
#
# You can remove the examples below

@article{Batut2018,
  doi = {10.1016/j.cels.2018.05.012},
  url = {https://doi.org/10.1016/j.cels.2018.05.012},
  year = {2018},
  month = jun,
  publisher = {Elsevier {BV}},
  volume = {6},
  number = {6},
  pages = {752--758.e1},
  author = {B{\\'{e}}r{\\'{e}}nice Batut and Saskia Hiltemann and Andrea Bagnacani and Dannon Baker and Vivek Bhardwaj and
           Clemens Blank and Anthony Bretaudeau and Loraine Brillet-Gu{\\'{e}}guen and Martin {\\v{C}}ech and John Chilton
           and Dave Clements and Olivia Doppelt-Azeroual and Anika Erxleben and Mallory Ann Freeberg and Simon Gladman and
           Youri Hoogstrate and Hans-Rudolf Hotz and Torsten Houwaart and Pratik Jagtap and Delphine Larivi{\\`{e}}re and
           Gildas Le Corguill{\\'{e}} and Thomas Manke and Fabien Mareuil and Fidel Ram{\\'{i}}rez and Devon Ryan and
           Florian Christoph Sigloch and Nicola Soranzo and Joachim Wolff and Pavankumar Videm and Markus Wolfien and
           Aisanjiang Wubuli and Dilmurat Yusuf and James Taylor and Rolf Backofen and Anton Nekrutenko and Bj\\"{o}rn Gr\\"{u}ning},
  title = {Community-Driven Data Analysis Training for Biology},
  journal = {Cell Systems}
}

@online{gtn-website,
  author = {GTN community},
  title = {GTN Training Materials: Collection of tutorials developed and maintained by the worldwide Galaxy community},
  url = {https://training.galaxyproject.org},
  urldate = {2021-03-24}
}

"""

TUTO_HAND_ON_BODY_TEMPLATE = """

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{{ '{%' }} raw {{ '%}' }} `{{ '{%' }} cite Batut2018 {{ '%}' }}`{{ '{%' }} endraw {{ '%}' }}.
This will be rendered like this: {{ '{%' }} cite Batut2018 {{ '%}' }}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.

<!-- This is a comment. -->

**Please follow our
[tutorial to learn how to fill the Markdown]({{ '{{' }} site.baseurl {{ '}}' }}/topics/contributing/tutorials/\
create-new-tutorial-content/tutorial.html)**

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ '{{' }} page.zenodo_link {{ '}}' }}) or from
>    the shared data library (`GTN - Material` -> `{{ '{{' }} page.topic_name {{ '}}' }}`
>     -> `{{ '{{' }} page.title {{ '}}' }}`):
>
>    ```
>    {{ z_file_links }}
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {{ '{%' }} snippet faqs/galaxy/datasets_import_via_link.md {{ '%}' }}
>
>    {{ '{%' }} snippet faqs/galaxy/datasets_import_from_data_library.md {{ '%}' }}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {{ '{%' }} snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" {{ '%}' }}
>
> 5. Add to each database a tag corresponding to ...
>
>    {{ '{%' }} snippet faqs/galaxy/datasets_add_tag.md {{ '%}' }}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:

{{ body }}

## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
"""


class Tutorial:
    """Class to describe a training tutorial."""

    def __init__(self, training, topic, name="new_tuto", title="The new tutorial", zenodo_link=""):
        """Init a tutorial instance."""
        self.training = training
        self.topic = topic
        self.name = name
        self.title = title
        self.zenodo_link = zenodo_link
        self.zenodo_file_links = []
        self.questions = []
        self.objectives = []
        self.time = ""
        self.key_points = []
        self.contributors = []
        self.body = ""
        self.init_wf_fp = None
        self.init_wf_id = None
        self.hands_on = True
        self.slides = False
        self.set_dir_name()
        self.init_data_lib()
        self.body = templates.render(
            HANDS_ON_TOOL_BOX_TEMPLATE,
            **{"tool_name": "My Tool", "inputlist": get_empty_input(), "paramlist": get_empty_param()},
        )

    def init_from_kwds(self, kwds):
        """Init a tutorial instance from a kwds dictionary."""
        self.name = kwds["tutorial_name"]
        self.title = kwds["tutorial_title"]
        self.zenodo_link = kwds["zenodo_link"] if kwds["zenodo_link"] else ""
        self.questions = [
            "Which biological questions are addressed by the tutorial?",
            "Which bioinformatics techniques are important to know for this type of data?",
        ]
        self.objectives = [
            "The learning objectives are the goals of the tutorial",
            "They will be informed by your audience and will communicate to them and to yourself what you should focus on during the course",
            "They are single sentences describing what a learner should be able to do once they have completed the tutorial",
            "You can use Bloom's Taxonomy to write effective learning objectives",
        ]
        self.time = "3H"
        self.key_points = ["The take-home messages", "They will appear at the end of the tutorial"]
        self.contributors = ["contributor1", "contributor2"]
        self.init_wf_fp = kwds["workflow"]
        self.init_wf_id = kwds["workflow_id"]
        self.hands_on = kwds["hands_on"]
        self.slides = kwds["slides"]
        self.set_dir_name()
        self.init_data_lib()

    def init_from_existing_tutorial(self, tuto_name):
        """Init a tutorial instance from an existing tutorial (data library and tutorial.md)."""
        self.name = tuto_name
        self.set_dir_name()

        if not self.exists():
            raise Exception("The tutorial %s does not exists. It should be created" % self.name)

        # get the metadata information of the tutorial (from the top of the tutorial.md)
        with open(self.tuto_fp) as tuto_f:
            tuto_content = tuto_f.read()
        regex = r"^---\n(?P<metadata>[\s\S]*)\n---(?P<body>[\s\S]*)"
        tuto_split_regex = re.search(regex, tuto_content)
        if not tuto_split_regex:
            raise Exception("No metadata found at the top of the tutorial")
        metadata = yaml.safe_load(tuto_split_regex.group("metadata"))
        self.title = metadata["title"]
        self.zenodo_link = metadata["zenodo_link"]
        self.questions = metadata["questions"]
        self.objectives = metadata["objectives"]
        self.time_estimation = metadata["time_estimation"]
        self.key_points = metadata["key_points"]
        self.contributors = metadata["contributors"]

        # the tutorial content
        self.body = tuto_split_regex.group("body")

        # get the data library
        self.init_data_lib()

    def init_data_lib(self):
        """Init the data library dictionary."""
        if os.path.exists(self.data_lib_fp):
            self.data_lib = load_yaml(self.data_lib_fp)
        else:
            self.data_lib = collections.OrderedDict()
        # set default information
        self.data_lib.setdefault("destination", collections.OrderedDict())
        self.data_lib["destination"]["type"] = "library"
        self.data_lib["destination"]["name"] = "GTN - Material"
        self.data_lib["destination"]["description"] = "Galaxy Training Network Material"
        self.data_lib["destination"][
            "synopsis"
        ] = "Galaxy Training Network Material. See https://training.galaxyproject.org"
        self.data_lib.setdefault("items", [])
        self.data_lib.pop("libraries", None)
        # get topic or create new one
        topic = collections.OrderedDict()
        for item in self.data_lib["items"]:
            if item["name"] == self.topic.title:
                topic = item
        if not topic:
            self.data_lib["items"].append(topic)
            topic["name"] = self.topic.title
            topic["description"] = self.topic.summary
            topic["items"] = []
        # get tutorial or create new one
        self.tuto_data_lib = collections.OrderedDict()
        for item in topic["items"]:
            if item["name"] == self.title:
                self.tuto_data_lib = item
        if not self.tuto_data_lib:
            topic["items"].append(self.tuto_data_lib)
            self.tuto_data_lib["name"] = self.title
            self.tuto_data_lib["items"] = []

    # GETTERS
    def get_tuto_metata(self):
        """Return the string corresponding to the tutorial metadata."""
        metadata = collections.OrderedDict()
        metadata["title"] = self.title
        metadata["zenodo_link"] = self.zenodo_link
        metadata["questions"] = self.questions
        metadata["objectives"] = self.objectives
        metadata["time_estimation"] = self.time
        metadata["key_points"] = self.key_points
        metadata["contributors"] = self.contributors
        return yaml.safe_dump(metadata, indent=2, default_flow_style=False, default_style="", explicit_start=False)

    # SETTERS
    def set_dir_name(self):
        """Set the path to dir and files of a tutorial."""
        self.dir = os.path.join(self.topic.dir, "tutorials", self.name)
        self.tuto_fp = os.path.join(self.dir, "tutorial.md")
        self.bib_fp = os.path.join(self.dir, "tutorial.bib")
        self.slide_fp = os.path.join(self.dir, "slides.html")
        self.data_lib_fp = os.path.join(self.dir, "data-library.yaml")
        self.wf_dir = os.path.join(self.dir, "workflows")
        self.wf_fp = os.path.join(self.wf_dir, "main_workflow.ga")
        self.faq_dir = os.path.join(self.dir, "faqs")
        self.tour_dir = os.path.join(self.dir, "tours")
        # remove empty workflow file if there
        empty_wf_filepath = os.path.join(self.wf_dir, "empty_workflow.ga")
        if os.path.exists(empty_wf_filepath):
            os.remove(empty_wf_filepath)

    # TEST METHODS
    def exists(self):
        """Test if the tutorial exists."""
        return os.path.isdir(self.dir)

    def has_workflow(self):
        """Test if a workflow is provided for the tutorial."""
        return self.init_wf_fp or self.init_wf_id

    # EXPORT METHODS
    def export_workflow_file(self):
        """Copy or extract workflow file and add it to the tutorial directory."""
        if not os.path.exists(self.wf_dir):
            os.makedirs(self.wf_dir)
        if not os.path.exists(os.path.join(self.wf_dir, "index.md")):
            with open(os.path.join(self.wf_dir, "index.md"), "w") as handle:
                handle.write("---\nlayout: workflow-list\n---\n")
        if self.init_wf_fp:
            shutil.copy(self.init_wf_fp, self.wf_fp)
        elif self.init_wf_id:
            gi = GalaxyInstance(self.training.galaxy_url, key=self.training.galaxy_api_key)
            gi.workflows.export_workflow_to_local_path(self.init_wf_id, self.wf_fp, use_default_filename=False)

    # OTHER METHODS
    def get_files_from_zenodo(self):
        """Extract a list of URLs and dictionary describing the files from the JSON output of the Zenodo API."""
        z_record, req_res = get_zenodo_record(self.zenodo_link)

        self.zenodo_file_links = []
        if "files" not in req_res:
            raise ValueError("No files in the Zenodo record")

        files = []
        for f in req_res["files"]:
            file_dict = {"url": "", "src": "url", "ext": "auto", "info": self.zenodo_link}
            if "links" not in f and "self" not in f["links"]:
                raise ValueError("No link for file %s" % f)
            file_dict["url"] = f["links"]["self"]
            self.zenodo_file_links.append(f["links"]["self"])
            files.append(file_dict)

        return (files, z_record)

    def prepare_data_library_from_zenodo(self):
        """Get the list of URLs of the files on Zenodo, fill the data library, save it into the file."""
        self.zenodo_file_links = []
        if self.zenodo_link != "":
            files, z_record = self.get_files_from_zenodo()
            if z_record:
                # get current data library and/or previous data library for the tutorial
                # remove the latest tag of any existing library
                # remove the any other existing library
                current_data_lib = collections.OrderedDict()
                previous_data_lib = collections.OrderedDict()
                for item in self.tuto_data_lib["items"]:
                    if item["name"] == "DOI: 10.5281/zenodo.%s" % z_record:
                        current_data_lib = item
                    elif item["description"] == "latest":
                        previous_data_lib = item
                        previous_data_lib["description"] = ""
                if not current_data_lib:
                    current_data_lib["name"] = "DOI: 10.5281/zenodo.%s" % z_record
                    current_data_lib["description"] = "latest"
                    current_data_lib["items"] = []
                current_data_lib["items"] = files

                self.tuto_data_lib["items"] = [current_data_lib]
                if previous_data_lib:
                    self.tuto_data_lib["items"].append(previous_data_lib)
        save_to_yaml(self.data_lib, self.data_lib_fp)

    def write_hands_on_tutorial(self, add_z_file_links=True):
        """Write the content of the hands-on tutorial in the corresponding file."""
        if add_z_file_links:
            self.body = templates.render(
                TUTO_HAND_ON_BODY_TEMPLATE,
                **{"z_file_links": "\n>    ".join(self.zenodo_file_links), "body": self.body},
            )
        # write in the tutorial file with the metadata on the top
        metadata = self.get_tuto_metata()
        with open(self.tuto_fp, "w") as md:
            md.write(templates.render(TUTO_HAND_ON_TEMPLATE, **{"metadata": metadata, "body": self.body}))

        # create the bibliography file
        self.write_bibliography()

    def write_bibliography(self):
        """Write the content of the bibliography file for the tutorial."""
        with open(self.bib_fp, "w") as bib:
            bib.write(templates.render(TUTO_BIBLIOGRAPHY_TEMPLATE, **{"body": self.body}))

    def create_hands_on_tutorial(self, ctx):
        """Create tutorial structure from the workflow file (if it is provided)."""
        # load workflow and get hands-on body from the workflow
        if self.init_wf_id:
            if not self.training.galaxy_url:
                raise ValueError("No Galaxy URL given")
            self.body = get_hands_on_boxes_from_running_galaxy(
                self.init_wf_id, self.training.galaxy_url, self.training.galaxy_api_key
            )
        elif self.init_wf_fp:
            self.body = get_hands_on_boxes_from_local_galaxy(self.training.kwds, self.init_wf_fp, ctx)
        # write tutorial body
        self.write_hands_on_tutorial()

    def create_tutorial(self, ctx):
        """Create the skeleton of a new tutorial."""
        # create tuto folder and empty files
        os.makedirs(self.dir)
        os.makedirs(self.tour_dir)
        os.makedirs(self.wf_dir)

        # extract the data library from Zenodo and the links for the tutorial
        if self.zenodo_link != "":
            info("Create the data library from Zenodo")
            self.prepare_data_library_from_zenodo()

        # create tutorial skeleton from workflow and copy workflow file
        if self.hands_on:
            info("Create tutorial skeleton from workflow (if it is provided)")
            self.create_hands_on_tutorial(ctx)
            self.export_workflow_file()

        # create slide skeleton
        if self.slides:
            with open(self.slide_fp, "w") as slide_f:
                slide_f.write(templates.render(TUTO_SLIDES_TEMPLATE, **{"metadata": self.get_tuto_metata()}))

        # create the FAQ page
        os.makedirs(self.faq_dir)
        if not os.path.exists(os.path.join(self.faq_dir, "index.md")):
            with open(os.path.join(self.faq_dir, "index.md"), "w") as handle:
                handle.write("---\nlayout: faq-page\n---\n")


def get_zenodo_record(zenodo_link):
    """Get the content of a Zenodo record."""
    # get the record in the Zenodo link
    if "doi" in zenodo_link:
        z_record = zenodo_link.split(".")[-1]
    else:
        z_record = zenodo_link.split("/")[-1]
    # get JSON corresponding to the record from Zenodo API
    req = "https://zenodo.org/api/records/%s" % (z_record)
    try:
        r = requests.get(req)
        r.raise_for_status()
        req_res = r.json()
    except Exception as e:
        error(f"The Zenodo link ({zenodo_link}) seems invalid: {e}")
        req_res = {"files": []}
        z_record = None
    return (z_record, req_res)


def get_wf_inputs(step_inp):
    """Get the inputs from a workflow step and format them into a hierarchical dictionary."""
    inputs = {}
    for inp_n, inp in step_inp.items():
        if "|" in inp_n:
            repeat_regex = r"(?P<prefix>[^\|]*)_(?P<nb>\d+)\|(?P<suffix>.+).+"
            repeat_search = re.search(repeat_regex, inp_n)
            hier_regex = r"(?P<prefix>[^\|]*)\|(?P<suffix>.+)"
            hier_regex = re.search(hier_regex, inp_n)
            if repeat_search and repeat_search.start(0) <= hier_regex.start(0):
                inputs.setdefault(repeat_search.group("prefix"), {})
                inputs[repeat_search.group("prefix")].setdefault(
                    repeat_search.group("nb"), get_wf_inputs({hier_regex.group("suffix"): inp})
                )
            else:
                inputs.setdefault(hier_regex.group("prefix"), {})
                inputs[hier_regex.group("prefix")].update(get_wf_inputs({hier_regex.group("suffix"): inp}))
        else:
            inputs.setdefault(inp_n, inp)
    return inputs


def get_wf_param_values(init_params, inp_connections):
    """Get the param values from a workflow step and format them into a hierarchical dictionary."""
    if not isinstance(init_params, str) or '": ' not in init_params:
        form_params = init_params
    else:
        form_params = json.loads(init_params)
    if isinstance(form_params, dict):
        if "__class__" in form_params and (
            form_params["__class__"] == "RuntimeValue" or form_params["__class__"] == "ConnectedValue"
        ):
            form_params = inp_connections
        else:
            for p in form_params:
                inp = inp_connections[p] if p in inp_connections else {}
                form_params[p] = get_wf_param_values(form_params[p], inp)
    elif isinstance(form_params, list):
        json_params = form_params
        form_params = []
        for i, p in enumerate(json_params):
            inp = inp_connections[str(i)] if str(i) in inp_connections else {}
            form_params.append(get_wf_param_values(p, inp))
    elif isinstance(form_params, str) and '"' in form_params:
        form_params = form_params.replace('"', "")
    return form_params


def format_wf_steps(wf, gi):
    """Get a string with the hands-on boxes describing the different steps of the worklow."""
    body = ""
    steps = wf["steps"]

    for s in range(len(steps)):
        wf_step = steps[str(s)]

        # get params in workflow
        wf_param_values = {}
        if wf_step["tool_state"] and wf_step["input_connections"]:
            wf_param_values = get_wf_param_values(wf_step["tool_state"], get_wf_inputs(wf_step["input_connections"]))
        if not wf_param_values:
            continue
        # get tool description
        try:
            tool_desc = gi.tools.show_tool(wf_step["tool_id"], io_details=True)
        except Exception:
            tool_desc = {"inputs": []}
        # get formatted param description
        paramlist = ""
        for inp in tool_desc["inputs"]:
            if inp["name"].startswith("__"):
                continue
            tool_inp = ToolInput(inp, wf_param_values, steps, 1, should_be_there=True)
            paramlist += tool_inp.get_formatted_desc()
        # format the hands-on box
        body += templates.render(
            HANDS_ON_TOOL_BOX_TEMPLATE,
            **{"tool_name": wf_step["name"], "tool_id": wf_step["tool_id"], "paramlist": paramlist},
        )
    return body


def get_hands_on_boxes_from_local_galaxy(kwds, wf_filepath, ctx):
    """Server local Galaxy and get the workflow dictionary."""
    assert is_galaxy_engine(**kwds)
    runnable = for_path(wf_filepath)
    tuto_body = ""
    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([runnable]) as config:
            info("Status of installed repositories: %s" % config.gi.toolshed.get_repositories())
            workflow_id = config.workflow_id(wf_filepath)
            wf = config.gi.workflows.export_workflow_dict(workflow_id)
            tuto_body = format_wf_steps(wf, config.gi)
    return tuto_body


def get_hands_on_boxes_from_running_galaxy(wf_id, galaxy_url, galaxy_api_key):
    """Get the workflow dictionary from a running Galaxy instance with the workflow installed on it."""
    gi = GalaxyInstance(galaxy_url, key=galaxy_api_key)
    wf = gi.workflows.export_workflow_dict(wf_id)
    tuto_body = format_wf_steps(wf, gi)
    return tuto_body
