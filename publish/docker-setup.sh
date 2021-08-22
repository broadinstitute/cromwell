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

apt-get update
apt-get install \
    adoptopenjdk-11-hotspot \
    -y --no-install-recommends

# sbt launcher non-deb package installation instructions adapted from
# - https://github.com/sbt/sbt/releases/tag/v1.4.9
# - https://github.com/broadinstitute/scala-baseimage/pull/4/files
curl --location --fail --silent --show-error "https://github.com/sbt/sbt/releases/download/v1.4.9/sbt-1.4.9.tgz" |
    tar zxf - -C /usr/share
update-alternatives --install /usr/bin/sbt sbt /usr/share/sbt/bin/sbt 1

# Update sbt launcher
sbt -Dsbt.supershell=false -Dsbt.rootdir=true sbtVersion
