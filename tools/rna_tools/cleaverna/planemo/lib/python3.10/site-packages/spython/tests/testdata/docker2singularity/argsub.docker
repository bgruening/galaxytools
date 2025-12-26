ARG CUDA_VERSION=11.1.1
ARG OS_VERSION=20.04

FROM nvidia/cuda:${CUDA_VERSION}-cudnn8-devel-ubuntu${OS_VERSION}

LABEL maintainer="Dong Wang"


ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ARG DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get upgrade -y &&\
apt-get install -y wget python3-pip

RUN python3 -m pip install --upgrade pip

WORKDIR /workspace
ADD ./requirements.txt /workspace
RUN python3 -m pip install -r /workspace/requirements.txt and &&\
rm /workspace/requirements.txt

CMD ["/bin/bash"]
