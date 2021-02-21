#!/usr/bin/env bash

# Installs required dependencies inside the docker image used for publishing cromwell

set -eou pipefail

apt-get update
apt-get install \
    apt-transport-https \
    curl \
    git \
    gnupg \
    wget \
    -y --no-install-recommends

# setup install for adoptopenjdk
# https://adoptopenjdk.net/installation.html#linux-pkg-deb
wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add -
echo "deb https://adoptopenjdk.jfrog.io/adoptopenjdk/deb $(
        grep UBUNTU_CODENAME /etc/os-release | cut -d = -f 2
    ) main" |
    tee /etc/apt/sources.list.d/adoptopenjdk.list

# Install jq 1.6 to ensure --rawfile is supported
curl -L https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 -o /usr/bin/jq
chmod +x /usr/bin/jq

# Install sbt via https://www.scala-sbt.org/1.0/docs/Installing-sbt-on-Linux.html
echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list
curl -sL "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x2EE0EA64E40A89B84B2DF73499E82A75642AC823" |
    apt-key add

apt-get update
apt-get install \
    adoptopenjdk-11-hotspot \
    sbt \
    -y --no-install-recommends

# Update sbt
sbt -Dsbt.rootdir=true sbtVersion
