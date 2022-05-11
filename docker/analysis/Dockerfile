FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
 && apt-get -y install \
    python3 \
    python3-pip \
 && rm -rf /var/cache/apt/*

COPY bin /opt
