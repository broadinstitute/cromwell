#! /bin/bash

#TODO: remove when done testing
export BRANCH=db_perf_scripts

set -x

# Make sure ip forwarding is enabled by default
echo "net.ipv4.ip_forward = 1" > /etc/sysctl.conf

# Install Make
apt-get install make

# Install docker-compose
curl -L https://github.com/docker/compose/releases/download/1.22.0/docker-compose-$(uname -s)-$(uname -m) -o /usr/local/bin/docker-compose
chmod +x /usr/local/bin/docker-compose

# Create the directory where everything is going
mkdir /app
cd /app

# Download the docker-compose script and other needed files
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/$BRANCH/scripts/perf/vm_scripts/docker_compose.yml .
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/$BRANCH/scripts/perf/vm_scripts/cromwell-dashboard.json .
mkdir cromwell
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/$BRANCH/scripts/perf/vm_scripts/cromwell/cromwell.conf cromwell/
mkdir mysql
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/$BRANCH/scripts/perf/vm_scripts/mysql/init_user.sql mysql/


cd /app

# Get custom attributes from instance metadata
export CROMWELL_VERSION_TAG=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/cromwell_version" -H "Metadata-Flavor: Google")
export CROMWELL_PROJECT=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/cromwell_project" -H "Metadata-Flavor: Google")
export CROMWELL_EXECUTION_ROOT=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/cromwell_bucket" -H "Metadata-Flavor: Google")
export CROMWELL_GRAFANA_HOST=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/cromwell_grafana_host" -H "Metadata-Flavor: Google")
export CROMWELL_GRAFANA_BUCKET=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/cromwell_grafana_bucket" -H "Metadata-Flavor: Google")

docker-compose up -d
