THIS TOOL WILL NOT RUN AT PRESENT.

The tool is 'work in progress' and needs at least the following sorting out:

1. Execution environment

Current the xchem_deep_score.py code can be run in the informaticsmatters/deep-app-ubuntu-1604:latest
container (see instructions at the top of the python file for doing so). The Galaxy execution environment needs
to define to run as this docker container.
Alternatively a conda environment could potentially be created but the dependencies are very complex and
some components need to be built from source.
Details for the dependencies are mostly described in the GitHub repo for the docker image:
https://github.com/InformaticsMatters/dls-deep/tree/ubuntu

2. GPU availability

The code must run in an environment with a GPU and with the CUDA drivers.
The docker image mentioned above has everything that is needed and will run on a GPU enabled environment
(a special version of Docker on the host machine is needed that supports GPUs).

Only the predictions need a GPU. The prior and latter steps run on CPU. Without a GPU you can specify the --mock
option which uses random numbers for the predicted scores.

3. Associated Python scripts.

The docker image contains additional python scripts (primarily /root/train/fragalysis_test_files/predict.py)
that are needed. If not running in a container these will need to be made available.