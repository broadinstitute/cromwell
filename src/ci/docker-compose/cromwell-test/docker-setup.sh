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
    mysql-client \
    postgresql-client \
    python-dev \
    software-properties-common \
    sudo \

# install docker
curl -fsSL https://get.docker.com -o get-docker.sh
sh get-docker.sh

cat <<CONF >/etc/init/docker-chown.conf
start on startup
task
exec chown root:docker /var/run/docker.sock
CONF

# install sbt
echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823
apt-get update
apt-get install -y sbt

# upgrade python dependencies
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py
pip install --upgrade --force-reinstall pyopenssl

useradd hoggett
usermod -aG docker hoggett
echo "hoggett ALL=NOPASSWD: ALL" >> /etc/sudoers
mkdir -p /home/hoggett
chown hoggett:hoggett /home/hoggett
