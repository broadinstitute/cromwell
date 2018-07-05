#!/usr/bin/env bash

set -e

apt-get update

apt-get install -y gnupg

echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823

apt-get install -y apt-utils
apt-get install -y apt-transport-https

apt-get update

apt-get install -y sudo
apt-get install -y sbt
apt-get install -y mysql-client
apt-get install -y build-essential
apt-get install -y curl
apt-get install -y python-dev

curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py

pip install --upgrade --force-reinstall pyopenssl
