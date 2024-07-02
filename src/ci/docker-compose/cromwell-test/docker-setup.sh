#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

# Don't prompt for inputs
# https://manpages.ubuntu.com/manpages/focal/en/man7/debconf.7.html#unattended%20package%20installation
export DEBIAN_FRONTEND=noninteractive

# install first round of packages without dependencies or required for those with dependencies
apt-get update
apt-get install -y \
    apt-utils \
    apt-transport-https \
    build-essential \
    ca-certificates \
    curl \
    gnupg \
    gnupg-agent \
    gnupg2 \
    jq \
    mysql-client \
    postgresql-client \
    python3-dev \
    software-properties-common \
    sudo \
    wget \

# setup install for Temurin
# https://adoptium.net/installation/linux/
mkdir -p /etc/apt/keyrings
wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc
echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^UBUNTU_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list

# setup install for gcloud
# https://cloud.google.com/sdk/docs/install#deb
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" |
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg |
    apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

# setup install for docker
# https://docs.docker.com/engine/install/ubuntu/
apt-get remove -y docker docker-engine docker.io containerd runc || true
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"

# install packages that required setup
apt-get update
apt-get install -y \
    temurin-11-jdk \
    containerd.io \
    docker-ce \
    docker-ce-cli \
    google-cloud-sdk \

# remove downloaded archive files
# https://manpages.ubuntu.com/manpages/focal/en/man8/apt-get.8.html#description
apt-get clean

# install sbt launcher
# non-deb package installation instructions adapted from
# - https://github.com/sbt/sbt/releases/tag/v1.4.9
# - https://github.com/broadinstitute/scala-baseimage/pull/4/files
curl -L --silent "https://github.com/sbt/sbt/releases/download/v1.10.0/sbt-1.10.0.tgz" |
    tar zxf - -C /usr/share
update-alternatives --install /usr/bin/sbt sbt /usr/share/sbt/bin/sbt 1

# install docker compose
# https://docs.docker.com/compose/install/
curl \
  --location --fail --silent --show-error \
  "https://github.com/docker/compose/releases/download/1.28.5/docker-compose-$(uname -s)-$(uname -m)" \
  -o /usr/local/bin/docker-compose
chmod +x /usr/local/bin/docker-compose

# set python as python3
# https://manpages.ubuntu.com/manpages/focal/en/man1/update-alternatives.1.html#commands
update-alternatives --install /usr/bin/python python /usr/bin/python3 1

# upgrade python dependencies
# https://pip.pypa.io/en/stable/installing/#installing-with-get-pip-py
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py
pip3 install --upgrade --force-reinstall pyopenssl
