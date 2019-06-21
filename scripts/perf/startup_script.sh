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

# Make sure ip forwarding is enabled by default so that docker doesn't loses connectivity
echo "net.ipv4.ip_forward = 1" > /etc/sysctl.conf

# Set up env variables
addVar CROMWELL_BRANCH=$(extract_metadata CROMWELL_PERF_SCRIPTS_BRANCH)
addVar TEST_CASE_DIRECTORY=$(extract_metadata TEST_CASE_DIRECTORY)
addVar CLOUD_SQL_DB_USER=$(extract_metadata CROMWELL_DB_USER)
addVar CLOUD_SQL_DB_PASSWORD=$(extract_metadata CROMWELL_DB_PASS)
addVar CLOUD_SQL_INSTANCE=$(extract_metadata CLOUD_SQL_INSTANCE)
addVar CROMWELL_DOCKER_IMAGE=$(extract_metadata CROMWELL_DOCKER_IMAGE)
addVar CROMWELL_PROJECT=$(extract_metadata CROMWELL_PROJECT)
addVar CROMWELL_EXECUTION_ROOT=$(extract_metadata CROMWELL_BUCKET)
addVar GCS_REPORT_BUCKET=$(extract_metadata GCS_REPORT_BUCKET)
addVar GCS_REPORT_PATH=$(extract_metadata GCS_REPORT_PATH)
addVar BUILD_ID=$(extract_metadata BUILD_TAG)
addVar CLEAN_UP=$(extract_metadata CLEAN_UP)

addVar REPORT_BUCKET=cromwell-perf-test-reporting

addVar CROMWELL_ROOT=/app
addVar PERF_ROOT=${CROMWELL_ROOT}/scripts/perf
addVar TEST_WORKFLOW_ROOT=${PERF_ROOT}/test_cases

addVar CROMWELL_UNDER_TEST="localhost"

# Clone cromwell to get the perf scripts. Use https to avoid ssh fingerprint prompt when the script runs
git clone -b ${CROMWELL_BRANCH} --depth 1 --single-branch https://github.com/broadinstitute/cromwell.git ${CROMWELL_ROOT}

source ${PERF_ROOT}/helper.inc.sh

if [ -n "${TEST_CASE_DIRECTORY}" ]
then
    # Must be a directory
    TEST_WORKFLOW_DIR=${TEST_WORKFLOW_ROOT}/${TEST_CASE_DIRECTORY}
    
    # If it contains a custom cromwell configuration, use that 
    if [ -f "${TEST_WORKFLOW_DIR}/cromwell.conf" ]
    then
        addVar CROMWELL_CONF_DIR=${TEST_WORKFLOW_DIR}
        # copy the default one next to it so we can include it
        cp ${PERF_ROOT}/vm_scripts/cromwell/cromwell.conf ${CROMWELL_CONF_DIR}/cromwell_default.conf
    else
        # Otherwise use the default one
        addVar CROMWELL_CONF_DIR=${PERF_ROOT}/vm_scripts/cromwell
    fi
    
    # If it contains a custom centaur configuration, use that 
    if [ -f "${TEST_WORKFLOW_DIR}/centaur.conf" ]
    then
        addVar CENTAUR_CONF_DIR=${TEST_WORKFLOW_DIR}
        # copy the default one next to it so we can include it
        cp ${PERF_ROOT}/vm_scripts/centaur/centaur.conf ${CENTAUR_CONF_DIR}/centaur_default.conf
    else
        # Otherwise use the default one
        addVar CENTAUR_CONF_DIR=${PERF_ROOT}/vm_scripts/centaur
    fi
    
    addVar CENTAUR_TEST_FILE=$(ls ${TEST_WORKFLOW_DIR}/*.test | head)
    sed -i "s/\$BRANCH/${CROMWELL_BRANCH}/" ${CENTAUR_TEST_FILE}
fi

# Start cromwell and cloud sql proxy
docker-compose -f ${PERF_ROOT}/vm_scripts/docker-compose.yml up -d

if [ -n "$CENTAUR_TEST_FILE" ]
then
    wait_for_cromwell
    run_test
else
    echo "No workflow provided, shutting down"
fi

shutdown
