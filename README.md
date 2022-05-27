[![Build Status](https://travis-ci.org/bgruening/galaxytools.svg?branch=master)](https://travis-ci.org/bgruening/galaxytools)
[![DOI](https://zenodo.org/badge/5466/bgruening/galaxytools.svg)](https://zenodo.org/badge/latestdoi/5466/bgruening/galaxytools)



Galaxy Tool wrappers
====================

This repository contains a variety of different tools that can be installed and used inside the Galaxy. Many tools are already included in the [Galaxy Tool Shed](https://toolshed.g2.bx.psu.edu/), others needs some more [love](https://github.com/bgruening/galaxytools/issues).


Wrapping tools for use in Galaxy is easy! If you want to start please see the [Galaxy wiki](https://wiki.galaxyproject.org/Admin/Tools/AddToolTutorial) or get in contact. If you want to contribute to this repository please see our [contributing guidelines](CONTRIBUTING.md) and our list of [issues](https://github.com/bgruening/galaxytools/issues).


Highlights
----------

  * [ChemicalToolBox](./chemicaltoolbox) includes a lot of cheminformatic tools into Galaxy
  * [RNA tools](./tools/rna_tools) are integrated as part of the [de.NBI](http://www.denbi.de) RNA Bioinformatics Center
  * Genome Annotation tools


Other repositories with high quality tools
------------------------------------------

 * [IUC repo](https://github.com/galaxyproject/tools-iuc)
 * [DevTeam repo](https://github.com/galaxyproject/tools-devteam)
 * [Peter Cock's blast repo](https://github.com/peterjc/galaxy_blast)
 * [Peter Cock's pico_galaxy repo](https://github.com/peterjc/pico_galaxy)
 * [Biopython repo](https://github.com/biopython/galaxy_packages)
 * [Galaxy Proteomics repo](https://github.com/galaxyproteomics/tools-galaxyp)

Running the tests
====================

First, install Planemo. Run the commands at the project root directory.
```
python3 -m venv planemo
. planemo/bin/activate
pip install planemo
```

To run the tests for a specific tool (e.g., pca) in a specific folder (e.g. sklearn)
```
planemo test --biocontainers tools/sklearn/pca.xml
```

To run the tests for all tools in a specific folder (e.g. sklearn)
```
planemo test --biocontainers tools/sklearn
```

To run the tests for all tools in all folders 
```
planemo test --biocontainers tools
```

3 steps to get your tool into Galaxy - A real-world example
==============================================================

In this blog post, we will explain how you can get your software tool into a Galaxy server and with this, 
exposed to thousands of researchers. For this purpose, we will follow David’s steps to add the very generic 
UNIX `diff` tool to Galaxy.

The first step to getting your software tool deployed into a Galaxy instance is to develop a Conda package for it.
Conda is the *de facto* standard in many different communities to deploy software easily and reproducibly.
The European Galaxy team is heavily involved in the [conda-forge](https://conda-forge.org) and 
[Bioconda](https://www.nature.com/articles/s41592-018-0046-7) projects and Galaxy does have built-in support for 
both channels. If your software tool is from the Biomedical domain, we recommend the Bioconda channel. Otherwise,
create a Conda package for conda-forge. Here, David has created the following 
[Pull Request (PR) against the conda-forge repo](https://github.com/conda-forge/staged-recipes/pull/11170):

**Step 1 - the Conda package**: [https://github.com/conda-forge/staged-recipes/pull/11170](https://github.com/conda-forge/staged-recipes/pull/11170)

After merging, a [diffutils repository](https://github.com/conda-forge/diffutils-feedstock) is
created and the [Conda package](https://anaconda.org/conda-forge/diffutils) is available usually within 30 min.

The second step is to create the Galaxy wrapper. A Galaxy wrapper is a formal description of all inputs,
outputs and parameters of your tool, so that Galaxy can generate a GUI out of it and later a command to send 
to the cluster. You will find a tutorial on how to create such a wrapper in the 
[planemo](https://planemo.readthedocs.io/en/latest/writing.html) documentation. The community has created a 
few best-practices for [Galaxy wrapper development](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html) 
and we recommend to follow them as this will ensure your tools are high-quality and can be deployed at
the big public Galaxy servers. David has created the following PR was created against a public repository that 
collects a variety of different tools.

**Step 2 - the Galaxy wrapper**: [https://github.com/bgruening/galaxytools/pull/966](https://github.com/bgruening/galaxytools/pull/966)

We recommend the submission of your tool to one of the bigger community projects like the ones listed below.
This has the advantage that you will most likely get a review and can improve your tool, but also get some
infrastructure for automated testing and ToolShed deployment for free.

Other repositories with Galaxy tools:
* [IUC repo](https://github.com/galaxyproject/tools-iuc)
* [Björn Grüning's repo](https://github.com/bgruening/galaxytools)
* Peter Cock's repos:
   * [blast repo](https://github.com/peterjc/galaxy_blast)
   * [pico repo](https://github.com/peterjc/pico_galaxy)
   * [mira repo](https://github.com/peterjc/galaxy_mira)
* [Metabolomics](https://github.com/workflow4metabolomics/tools-metabolomics)
* [Proteomics](https://github.com/galaxyproteomics/tools-galaxyp)
* [Colibread Galaxy Tools](https://github.com/genouest/tools-colibread)
* [Greg von Kuster's repo](https://github.com/gregvonkuster/galaxy-csg)
* [EI repo](https://github.com/TGAC/earlham-galaxytools)
* [AAFC-MBB Canada repo](https://github.com/AAFC-MBB/Galaxy/tree/master/wrappers)
* [Mark Einon's repo](https://gitlab.com/einonm/galaxy-tools)
* [National Microbiology Laboratory's repo](https://github.com/phac-nml/galaxy_tools)

Once David’s Galaxy wrapper PR passed all tests and was merged, it was automatically pushed to the 
[Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/view/bgruening/diff/), an app store for Galaxy. 
From there, every Galaxy instance can install tools (apps).

Furthermore, a bot is automatically creating [(Bio)Containers](https://biocontainers.pro) (Docker, rkt and Singularity) by
tracking all Galaxy tools to ensure that a container exists for each tool. You can see the bot in action in the following PR:

Automatic containers: [https://github.com/BioContainers/multi-package-containers/pull/1236](https://github.com/BioContainers/multi-package-containers/pull/1236)

Last but not least, David wanted to get the Galaxy `diff` tool into the European Galaxy server.
For that, a new PR was created against the tool repository from usegalaxy-eu.

**Step 3 - request for installation**: [https://github.com/usegalaxy-eu/usegalaxy-eu-tools/pull/318](https://github.com/usegalaxy-eu/usegalaxy-eu-tools/pull/318)

Once this is merged, another bot installs all the new tools but also tool updates automatically every Saturday.
As a result, the installed `diff` tool can be used on the European Galaxy server following this link: https://usegalaxy.eu/root?tool_id=diff

That's it - 3 steps to get your tool exposed to thousands of researchers!
