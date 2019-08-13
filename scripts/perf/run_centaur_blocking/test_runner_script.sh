#! /bin/bash

# Deliberately don't set these options. The script must continue to run (and clean itself up) even if the centaur tests
# fail.
#set -euo pipefail
set +x

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

# Extract env variables from instance metadata:
addVar TEST_CASE_DIRECTORY=$(extract_metadata TEST_CASE_DIRECTORY)
addVar CROMWELL_UNDER_TEST=$(extract_metadata CROMWELL_UNDER_TEST)
addVar CROMWELL_PROJECT=$(extract_metadata CROMWELL_PROJECT)
addVar GCS_REPORT_BUCKET=$(extract_metadata GCS_REPORT_BUCKET)
addVar GCS_REPORT_PATH=$(extract_metadata GCS_REPORT_PATH)
addVar BUILD_ID=$(extract_metadata BUILD_TAG)
addVar CROMWELL_PERF_SCRIPTS_BRANCH=$(extract_metadata CROMWELL_PERF_SCRIPTS_BRANCH)

addVar CROMWELL_ROOT=/app
addVar PERF_ROOT=${CROMWELL_ROOT}/scripts/perf

# Clone cromwell to get the perf scripts. Use https to avoid ssh fingerprint prompt when the script runs
git clone -b ${CROMWELL_PERF_SCRIPTS_BRANCH} --depth 1 --single-branch https://github.com/broadinstitute/cromwell.git ${CROMWELL_ROOT}

source ${PERF_ROOT}/helper.inc.sh

addVar CLEAN_UP=true

addVar REPORT_BUCKET="cromwell-perf-test-reporting"
addVar TEST_WORKFLOW_ROOT="${PERF_ROOT}/test_cases"
addVar REPORT_URL="$(strip_trailing_slash "gs://${GCS_REPORT_BUCKET}/${GCS_REPORT_PATH}")"

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
    sed -i "s/\$BRANCH/${CROMWELL_PERF_SCRIPTS_BRANCH}/" ${CENTAUR_TEST_FILE}
fi

if [ -n "${CENTAUR_TEST_FILE}" ]
then
    wait_for_cromwell
    run_test
else
    echo "No workflow provided, shutting down"
fi

export_centaur_logs
self_destruct_instance
