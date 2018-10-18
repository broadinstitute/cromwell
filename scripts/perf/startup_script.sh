#! /bin/bash

set -x

### /!\ This script assumes docker and docker compose are already installed on the host

# Utility function to extract values from instance metadata
extract_metadata() {
  curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1" -H "Metadata-Flavor: Google"
}

# Make sure ip forwarding is enabled by default so that docker doesn't loses connectivity
echo "net.ipv4.ip_forward = 1" > /etc/sysctl.conf

# Set up env variables
export CROMWELL_BRANCH=$(extract_metadata CROMWELL_BRANCH_NAME)
export WORKFLOW_SCRIPT_NAME=$(extract_metadata WORKFLOW_SCRIPT)
export BUILD_ID=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
export CLOUD_SQL_DB_USER=$(extract_metadata CROMWELL_DB_USER)
export CLOUD_SQL_DB_PASSWORD=$(extract_metadata CROMWELL_DB_PASS)
export CLOUD_SQL_INSTANCES=$(extract_metadata CLOUD_SQL_INSTANCE)
export CROMWELL_IMAGE=$(extract_metadata CROMWELL_IMAGE)
export CROMWELL_PROJECT=$(extract_metadata CROMWELL_PROJECT)
export CROMWELL_EXECUTION_ROOT=$(extract_metadata CROMWELL_BUCKET)
export CROMWELL_STATSD_HOST=$(extract_metadata CROMWELL_STATSD_HOST)
export CROMWELL_STATSD_PORT=$(extract_metadata CROMWELL_STATSD_PORT)
export CROMWELL_STATSD_PREFIX=$(extract_metadata CROMWELL_STATSD_PREFIX)
export DO_SHUTDOWN=$(extract_metadata DO_SHUTDOWN)
# Use the instance name as statsd prefix to avoid metrics collisions
export CROMWELL_STATSD_PREFIX=${BUILD_ID}
export REPORT_BUCKET=cromwell-perf-test-reporting

export CROMWELL_ROOT=/app
export PERF_ROOT=${CROMWELL_ROOT}/scripts/perf
export TEST_WORKFLOW_ROOT=${PERF_ROOT}/test_cases

export CENTAUR_TEST_NAME=unknown-test
export CROMWELL_VERSION=unknown-version
export 

# Clone cromwell to get the perf scripts. Use https to avoid ssh fingerprint prompt when the script runs
git clone -b ${CROMWELL_BRANCH} --depth 1 --single-branch https://github.com/broadinstitute/cromwell.git ${CROMWELL_ROOT}

source ${PERF_ROOT}/helper.inc.sh

if [ -n "$WORKFLOW_SCRIPT_NAME" ]
then
    # Must be a directory
    TEST_WORKFLOW_DIR=${TEST_WORKFLOW_ROOT}/${WORKFLOW_SCRIPT_NAME}
    
    # If it contains a custom cromwell configuration, use that 
    if [ -f "${TEST_WORKFLOW_DIR}/cromwell.conf" ]
    then
        export CROMWELL_CONF_DIR=${TEST_WORKFLOW_DIR}
        # copy the default one next to it so we can include it
        cp ${PERF_ROOT}/vm_scripts/cromwell/cromwell.conf ${CROMWELL_CONF_DIR}/cromwell_default.conf
    else
        # Otherwise use the default one
        export CROMWELL_CONF_DIR=${PERF_ROOT}/vm_scripts/cromwell
    fi
    
    # If it contains a custom centaur configuration, use that 
    if [ -f "${TEST_WORKFLOW_DIR}/centaur.conf" ]
    then
        export CENTAUR_CONF_DIR=${TEST_WORKFLOW_DIR}
        # copy the default one next to it so we can include it
        cp ${PERF_ROOT}/vm_scripts/centaur/centaur.conf ${CENTAUR_CONF_DIR}/centaur_default.conf
    else
        # Otherwise use the default one
        export CENTAUR_CONF_DIR=${PERF_ROOT}/vm_scripts/centaur
    fi
    
    export CENTAUR_TEST_FILE=$(ls ${TEST_WORKFLOW_DIR}/*.test | head)
    sed -i "s/\$BRANCH/${CROMWELL_BRANCH}/" ${CENTAUR_TEST_FILE}
fi

# Clone the CloudSQL DB
gcloud --project broad-dsde-cromwell-perf sql instances clone cromwell-perf-testing-base-09-24-18 cromwell-db-${BUILD_ID}
gcloud --project broad-dsde-cromwell-perf sql users create cromwell --instance=cromwell-db-${BUILD_ID} --password=${CLOUD_SQL_DB_PASSWORD}

set_up

# Start cromwell and cloud sql proxy
prepare_statsd_proxy
docker-compose -f ${PERF_ROOT}/vm_scripts/docker-compose.yml up -d

if [ -n "$CENTAUR_TEST_FILE" ]
then
    wait_for_cromwell
    run_test
else
    echo "No workflow provided, shutting down"
fi

shutdown
