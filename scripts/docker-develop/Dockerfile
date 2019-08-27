FROM openjdk:8u171

# docker build -t vanessa/cromwell-dev .

# Development environment for Cromwell that includes:
#
#   Scala 2.12
#   SBT 1.x
#   Java 8
#   Git

# Env variables
ENV SCALA_VERSION 2.12.9
ENV SBT_VERSION 1.2.8

# Scala expects this file
RUN touch /usr/lib/jvm/java-8-openjdk-amd64/release

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
  sbt sbtVersion

# Instruct user to add code here during development
RUN useradd pigman && \
    mkdir -p /code

WORKDIR /code
