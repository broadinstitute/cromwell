#!/usr/bin/env bash

# Installs required dependencies inside the docker image used for publishing cromwell

set -eo pipefail # SDKMAN relies on testing/setting unbound variables so no `u`

apt-get update
apt-get install \
    apt-transport-https \
    curl \
    git \
    gnupg \
    wget \
    ca-certificates \
    unzip \
    zip \
    -y --no-install-recommends

# Recommended by SBT
# https://www.scala-sbt.org/1.x/docs/Installing-sbt-on-Linux.html#Installing+from+SDKMAN
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk version
sdk install java $(sdk list java | grep -o "\b11\.[0-9]*\.[0-9]*\-tem" | head -1) # latest `11` build of Temurin
java -version
sdk install sbt

# Install jq 1.6 to ensure --rawfile is supported
curl -L https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 -o /usr/bin/jq
chmod +x /usr/bin/jq

# Update sbt launcher
sbt -Dsbt.supershell=false -Dsbt.rootdir=true sbtVersion
