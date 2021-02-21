FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -q && \
    apt-get upgrade -qq && \
    apt-get install -y curl wget gnupg gnupg2 && \
    rm -rf /var/lib/apt/lists/*

# docker build -t vanessa/cromwell-dev .

# Development environment for Cromwell that includes:
#
#   Scala 2.12
#   SBT 1.x
#   Java 11
#   Git

# Env variables
ENV SCALA_VERSION 2.12.12
ENV SBT_VERSION 1.4.7

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
    curl -fsL https://downloads.typesafe.com/scala/$SCALA_VERSION/scala-$SCALA_VERSION.tgz | tar xfz - -C /opt/ && \
    echo >> /home/pigman/.bashrc && \
    echo "export PATH=/opt/scala-$SCALA_VERSION/bin:$PATH" >> /home/pigman/.bashrc

#
## sbt
#

RUN curl -L -o sbt-$SBT_VERSION.deb https://dl.bintray.com/sbt/debian/sbt-$SBT_VERSION.deb && \
  dpkg -i sbt-$SBT_VERSION.deb && \
  rm sbt-$SBT_VERSION.deb && \
  apt-get update && \
  apt-get install sbt && \
  sbt -Dsbt.rootdir=true sbtVersion

# Instruct user to add code here during development
RUN useradd pigman && \
    mkdir -p /code

WORKDIR /code
