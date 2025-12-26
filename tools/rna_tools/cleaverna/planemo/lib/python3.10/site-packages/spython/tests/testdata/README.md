# Test Data

This folder contains test data for Singularity Python.

 - [Singularity](Singularity) and [Docker](Docker) are generic recipes used to run the tests in [one folder up](../). They are not inended to be converted between one another.
 - [singularity2docker](singularity2docker) is a folder of docker recipes (`*.docker`) and Singularity recipes (`*.def`) that are tested for conversion *from* Singularity to Docker.
 - [docker2singularity](docker2singularity) is a folder of Singularity recipes (`*.def`) and docker recipes (`*.docker`) that are tested for conversion *from* Docker to Singularity.

To add a new pair of recipes to either folder, simply write a .def and .docker file with the same name. They will be tested by [test_conversion.py](../test_conversion.py).
