FROM ubuntu:16.04

RUN apt-get update && \
  apt-get install -y \
    python \
    python-pip \
    git && \
  rm -rf /var/lib/apt/lists/*

RUN pip install \
    biopython==1.69

RUN cd / && git clone https://github.com/shu65/astral_roc.git 

ENV PYTHONPATH=$PYTHONPATH:/astral_roc
WORKDIR /astral_roc/src/