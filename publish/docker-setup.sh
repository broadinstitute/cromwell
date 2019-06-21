#!/usr/bin/env bash

# Installs required dependencies inside the docker image used for publishing cromwell

set -eou pipefail

apt-get update
apt-get install \
    apt-transport-https \
    curl \
    git \
    gnupg \
    openjdk-8-jdk \
    -y --no-install-recommends

# Install jq 1.6 to ensure --rawfile is supported
curl -L https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 -o /usr/bin/jq
chmod +x /usr/bin/jq

# Install sbt via https://www.scala-sbt.org/1.0/docs/Installing-sbt-on-Linux.html
echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823
apt-get update
apt-get install sbt -y --no-install-recommends

# Update sbt
sbt sbtVersion
