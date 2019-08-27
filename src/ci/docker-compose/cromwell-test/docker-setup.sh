#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

apt-get update

# install mysql and other dependencies
apt-get install -y \
    apt-utils \
    apt-transport-https \
    build-essential \
    ca-certificates \
    curl \
    gnupg \
    gnupg2 \
    jq \
    mysql-client \
    postgresql-client \
    python-dev \
    software-properties-common \
    sudo \

# install sbt
echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823
apt-get update
apt-get install -y sbt

# upgrade python dependencies
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py
pip install --upgrade --force-reinstall pyopenssl

# install gcloud
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
apt-get update
apt-get install -y google-cloud-sdk

useradd hoggett
echo "hoggett ALL=NOPASSWD: ALL" >> /etc/sudoers
mkdir -p /home/hoggett
chown hoggett:hoggett /home/hoggett
