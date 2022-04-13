FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -q && \
    apt-get upgrade -qq && \
    apt-get install -y curl wget gnupg gnupg2 && \
    rm -rf /var/lib/apt/lists/*

# docker build -t vanessa/cromwell-dev .

# Development environment for Cromwell that includes:
#
#   Scala 2.13
#   SBT 1.x
#   Java 11
#   Git

# Env variables
ENV SCALA_VERSION 2.13.8
ENV SBT_VERSION 1.5.5

#
## AdoptOpenJDK Hotspot
#
# https://adoptopenjdk.net/installation.html#linux-pkg-deb
RUN \
wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add - && \
echo \
"deb https://adoptopenjdk.jfrog.io/adoptopenjdk/deb $(grep UBUNTU_CODENAME /etc/os-release | cut -d = -f 2) main" | \
tee /etc/apt/sources.list.d/adoptopenjdk.list && \
apt-get update && \
apt-get install -y adoptopenjdk-11-hotspot

#
## Scala
#

RUN mkdir -p /home/pigman && \
    curl \
    --location --fail --silent --show-error \
    https://downloads.typesafe.com/scala/$SCALA_VERSION/scala-$SCALA_VERSION.tgz | tar xfz - -C /opt/ && \
    echo >> /home/pigman/.bashrc && \
    echo "export PATH=/opt/scala-$SCALA_VERSION/bin:$PATH" >> /home/pigman/.bashrc

#
## sbt
#

# non-deb package installation instructions adapted from
# - https://github.com/sbt/sbt/releases/tag/v1.4.9
# - https://github.com/broadinstitute/scala-baseimage/pull/4/files
RUN curl \
    --location --fail --silent --show-error \
    "https://github.com/sbt/sbt/releases/download/v$SBT_VERSION/sbt-$SBT_VERSION.tgz" | \
    tar zxf - -C /usr/share && \
    update-alternatives --install /usr/bin/sbt sbt /usr/share/sbt/bin/sbt 1 && \
    sbt -Dsbt.supershell=false -Dsbt.rootdir=true sbtVersion

# Instruct user to add code here during development
RUN useradd pigman && \
    mkdir -p /code

WORKDIR /code
