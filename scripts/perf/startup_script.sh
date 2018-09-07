#! /bin/bash

### /!\ This script assumes docker and docker compose are already installed on the host

# Utility function to extract values from instance metadata
extract_metadata() {
  curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1" -H "Metadata-Flavor: Google"
}

# Get Build ID
export BUILD_ID=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")


# Get user/password
export CLOUD_SQL_DB_USER=$(extract_metadata CROMWELL_DB_USER)
export CLOUD_SQL_DB_PASSWORD=$(extract_metadata CROMWELL_DB_PASS)

gcloud --project broad-dsde-cromwell-perf sql instances clone cromwell-perf-testing-base cromwell-db-$BUILD_ID
gcloud --project broad-dsde-cromwell-perf sql users create cromwell --instance=cromwell-db-$BUILD_ID --password=$CLOUD_SQL_DB_PASSWORD

set -x

# Make sure ip forwarding is enabled by default
echo "net.ipv4.ip_forward = 1" > /etc/sysctl.conf

# Create the directory where everything is going
mkdir /app
cd /app

# Extract branch name
export CROMWELL_BRANCH=$(extract_metadata CROMWELL_BRANCH_NAME)

# Download the docker-compose script and cromwell configuration
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/${CROMWELL_BRANCH}/scripts/perf/vm_scripts/docker-compose.yml -o docker-compose.yml
mkdir cromwell
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/${CROMWELL_BRANCH}/scripts/perf/vm_scripts/cromwell/cromwell.conf -o cromwell/cromwell.conf

# Get custom attributes from instance metadata
export CLOUD_SQL_INSTANCES=$(extract_metadata CLOUD_SQL_INSTANCE)
export CROMWELL_VERSION_TAG=$(extract_metadata CROMWELL_VERSION)
export CROMWELL_PROJECT=$(extract_metadata CROMWELL_PROJECT)
export CROMWELL_EXECUTION_ROOT=$(extract_metadata CROMWELL_BUCKET)
export CROMWELL_STATSD_HOST=$(extract_metadata CROMWELL_STATSD_HOST)
export CROMWELL_STATSD_PORT=$(extract_metadata CROMWELL_STATSD_PORT)
export CROMWELL_STATSD_PREFIX=$(extract_metadata CROMWELL_STATSD_PREFIX)

# Use the instance name as statsd prefix to avoid metrics collisions
export CROMWELL_STATSD_PREFIX=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")

docker-compose up -d

echo "Trying to ping Cromwell"

until [[ $(curl -X GET "http://localhost:8000/engine/v1/version" -H "accept: application/json") = *"cromwell"* ]]; do
	echo "No connection to Cromwell yet!"
	sleep 5
done

echo "Connection to Cromwell success!"

mkdir /workflow_files
cd /workflow_files

curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/develop/centaur/src/main/resources/standardTestCases/hello/hello.wdl -o workflow.wdl
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/develop/centaur/src/main/resources/standardTestCases/hello/hello.inputs -o workflow_inputs.json

curl -X POST "http://localhost:8000/api/workflows/v1" -H "accept: application/json" -H "Content-Type: multipart/form-data" -F "workflowSource=@workflow.wdl" -F "workflowInputs=@workflow_inputs.json;type=application/json"
