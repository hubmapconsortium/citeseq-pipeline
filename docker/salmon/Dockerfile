FROM ubuntu:22.04

WORKDIR /opt

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
 && apt-get -y install \
    python3 \
    python3-pip \
    salmon \
    git-all \
 && rm -rf /var/cache/apt/*

COPY docker/salmon/requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN mkdir /out

COPY bin /opt
RUN chmod -R +rwx /opt
