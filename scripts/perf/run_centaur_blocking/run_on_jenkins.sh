#!/usr/bin/env bash

set -euo pipefail

source scripts/perf/helper.inc.sh

export CROMWELL_ROOT=$(pwd)
export PERF_ROOT=${CROMWELL_ROOT}/scripts/perf
export TEST_WORKFLOW_ROOT=${PERF_ROOT}/test_cases

if [ -n "${TEST_CASE_DIRECTORY}" ]
then
    # Must be a directory
    TEST_WORKFLOW_DIR=${TEST_WORKFLOW_ROOT}/${TEST_CASE_DIRECTORY}

    # If it contains a custom cromwell configuration, use that
    if [ -f "${TEST_WORKFLOW_DIR}/cromwell.conf" ]
    then
        export CROMWELL_CONF_DIR="${TEST_WORKFLOW_DIR}"
        # copy the default one next to it so we can include it
        cp ${PERF_ROOT}/vm_scripts/cromwell/cromwell.conf ${CROMWELL_CONF_DIR}/cromwell_default.conf
    else
        # Otherwise use the default one
        export CROMWELL_CONF_DIR="${PERF_ROOT}/vm_scripts/cromwell"
    fi

    # If it contains a custom centaur configuration, use that
    if [ -f "${TEST_WORKFLOW_DIR}/centaur.conf" ]
    then
        export CENTAUR_CONF_DIR="${TEST_WORKFLOW_DIR}"
        # copy the default one next to it so we can include it
        cp ${PERF_ROOT}/vm_scripts/centaur/centaur.conf ${CENTAUR_CONF_DIR}/centaur_default.conf
    else
        # Otherwise use the default one
        export CENTAUR_CONF_DIR="${PERF_ROOT}/vm_scripts/centaur"
    fi

    export CENTAUR_TEST_FILE=$(ls ${TEST_WORKFLOW_DIR}/*.test | head)
    sed -i "s/\$BRANCH/${CROMWELL_BRANCH}/" "${CENTAUR_TEST_FILE}"
fi

if [ -n "${CENTAUR_TEST_FILE}" ]
then
    custom_wait_for_cromwell
    run_test
else
    echo "No workflow provided, shutting down"
fi
