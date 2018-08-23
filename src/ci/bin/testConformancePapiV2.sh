#!/usr/bin/env bash

set -e
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_conformance_environment

cromwell::build::assemble_jars

CENTAUR_CWL_RUNNER_MODE="papi"
GOOGLE_AUTH_MODE="service-account"
GOOGLE_SERVICE_ACCOUNT_JSON="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-centaur-service-account.json"
GOOGLE_SERVICE_ACCOUNT_JSON_REQUESTER_PAYS="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-centaur-requester-pays-service-account.json"
PAPI_INPUT_GCS_PREFIX=gs://centaur-cwl-conformance/cwl-inputs/

# Export variables used in conf files and commands
export CENTAUR_CWL_RUNNER_MODE
export GOOGLE_AUTH_MODE
export GOOGLE_SERVICE_ACCOUNT_JSON
export GOOGLE_SERVICE_ACCOUNT_JSON_REQUESTER_PAYS
export PAPI_INPUT_GCS_PREFIX

shutdown_cromwell() {
    if [ -n "${CROMWELL_PID+set}" ]; then
        cromwell::build::kill_tree "${CROMWELL_PID}"
    fi
}

cromwell::build::add_exit_function shutdown_cromwell

# Start the Cromwell server in the directory containing input files so it can access them via their relative path
cd "${CROMWELL_BUILD_CWL_TEST_RESOURCES}"

# Turn off call caching as hashing doesn't work since it sees local and not GCS paths.
# CWL conformance uses alpine images that do not have bash.
java \
    -Xmx2g \
    -Dconfig.file="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2_application.conf" \
    -Dcall-caching.enabled=false \
    -Dsystem.job-shell=/bin/sh \
    -jar "${CROMWELL_BUILD_JAR}" \
    server &

CROMWELL_PID=$!

sleep 30

cat <<JSON >"${CROMWELL_BUILD_CWL_TEST_INPUTS}"
{
    "cwl_conformance_test.cwl_dir": "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}",
    "cwl_conformance_test.test_result_output": "${CROMWELL_BUILD_CWL_TEST_OUTPUT}",
    "cwl_conformance_test.centaur_cwl_runner": "${CROMWELL_BUILD_CWL_TEST_RUNNER}",
    "cwl_conformance_test.conformance_expected_failures":
        "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_conformance_expected_failures.txt",
    "cwl_conformance_test.timeout": 1200
}
JSON

java \
    -Xmx2g \
    -Dbackend.providers.Local.config.concurrent-job-limit="${CROMWELL_BUILD_CWL_TEST_PARALLELISM}" \
    -jar "${CROMWELL_BUILD_JAR}" \
    run "${CROMWELL_BUILD_CWL_TEST_WDL}" \
    -i "${CROMWELL_BUILD_CWL_TEST_INPUTS}"

cd "${CROMWELL_BUILD_ROOT_DIRECTORY}"

cromwell::build::generate_code_coverage
