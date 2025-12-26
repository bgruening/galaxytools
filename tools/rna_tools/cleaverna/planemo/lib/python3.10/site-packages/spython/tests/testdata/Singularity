Bootstrap: docker
From: continuumio/miniconda3

%runscript
    exec /opt/conda/bin/spython "$@"

%labels
    maintainer vsochat@stanford.edu

%post
    apt-get update && apt-get install -y git

# Dependencies
    cd /opt
    git clone https://www.github.com/singularityhub/singularity-cli
    cd singularity-cli
    /opt/conda/bin/pip install setuptools
    /opt/conda/bin/python setup.py install
