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

# setup install for adoptopenjdk
# https://adoptopenjdk.net/installation.html#linux-pkg-deb
wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add -
echo "deb https://adoptopenjdk.jfrog.io/adoptopenjdk/deb $(
        grep UBUNTU_CODENAME /etc/os-release | cut -d = -f 2
    ) main" |
    tee /etc/apt/sources.list.d/adoptopenjdk.list

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
    adoptopenjdk-11-hotspot \
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
curl -L --silent "https://github.com/sbt/sbt/releases/download/v1.4.9/sbt-1.4.9.tgz" |
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

# create a non-root user with access to sudo (but not `sudo -u` / `sudo -g`)
# https://www.sudo.ws/man/1.9.4/sudoers.man.html#Runas_Spec
# https://www.sudo.ws/man/1.9.4/sudoers.man.html#NOPASSWD
# https://www.sudo.ws/man/1.9.4/sudoers.man.html#EXAMPLES
useradd hoggett
echo "hoggett ALL=NOPASSWD: ALL" >> /etc/sudoers
mkdir -p /home/hoggett
chown hoggett:hoggett /home/hoggett
