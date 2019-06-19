#! /bin/bash

set -x

### /!\ This script assumes docker and docker compose are already installed on the host

# Utility function to extract values from instance metadata
extract_metadata() {
  curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1" -H "Metadata-Flavor: Google"
}

# Exports an env var and also adds it to the root bashrc. This way if there is a need to ssh onto the machine
# for debugging one will have the env variables already set when using root
addVar() {
  export $1
  echo "export $1" >> /root/.bashrc
}

# Make sure ip forwarding is enabled by default so that docker doesn't lose connectivity
echo "net.ipv4.ip_forward = 1" > /etc/sysctl.conf

# Set up env variables
addVar CROMWELL_BRANCH=$(extract_metadata CROMWELL_PERF_SCRIPTS_BRANCH)
addVar CLOUD_SQL_DB_USER=$(extract_metadata CROMWELL_DB_USER)
addVar CLOUD_SQL_DB_PASSWORD=$(extract_metadata CROMWELL_DB_PASS)
addVar CLOUD_SQL_INSTANCE=$(extract_metadata CLOUD_SQL_INSTANCE)
addVar CROMWELL_DOCKER_IMAGE=$(extract_metadata CROMWELL_DOCKER_IMAGE)
addVar CROMWELL_PROJECT=$(extract_metadata CROMWELL_PROJECT)
addVar CROMWELL_EXECUTION_ROOT=$(extract_metadata CROMWELL_BUCKET)
addVar CROMWELL_STATSD_HOST=$(extract_metadata CROMWELL_STATSD_HOST)
addVar CROMWELL_STATSD_PORT=$(extract_metadata CROMWELL_STATSD_PORT)
addVar BUILD_ID=$(extract_metadata BUILD_TAG)

# Use the instance name as statsd prefix to avoid metrics collisions
addVar CROMWELL_STATSD_PREFIX=${BUILD_ID}
addVar REPORT_BUCKET=cromwell-perf-test-reporting

addVar CROMWELL_ROOT=/app
addVar PERF_ROOT=${CROMWELL_ROOT}/scripts/perf
addVar TEST_WORKFLOW_ROOT=${PERF_ROOT}/test_cases

# Clone cromwell to get the perf scripts. Use https to avoid ssh fingerprint prompt when the script runs
git clone -b ${CROMWELL_BRANCH} --depth 1 --single-branch https://github.com/broadinstitute/cromwell.git ${CROMWELL_ROOT}

source ${PERF_ROOT}/helper.inc.sh

addVar CROMWELL_CONF_DIR=${PERF_ROOT}/vm_scripts/cromwell

# Start cromwell and cloud sql proxy
prepare_statsd_proxy
docker-compose -f ${PERF_ROOT}/vm_scripts/docker-compose.yml up -d
