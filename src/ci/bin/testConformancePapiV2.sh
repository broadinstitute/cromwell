#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_conformance_environment

# Override of the default sbt assembly command which is just `assembly`.
# The conformance runs only need these two subprojects so save a couple of minutes and skip the rest.
CROMWELL_SBT_ASSEMBLY_COMMAND="server/assembly centaurCwlRunner/assembly"
cromwell::build::assemble_jars

CENTAUR_CWL_RUNNER_MODE="papi"
GOOGLE_AUTH_MODE="service-account"
PAPI_INPUT_GCS_PREFIX=gs://centaur-cwl-conformance/cwl-inputs/

# Export variables used in conf files and commands
export CENTAUR_CWL_RUNNER_MODE
export GOOGLE_AUTH_MODE
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
    -Dconfig.file="${CROMWELL_BUILD_CROMWELL_CONFIG}" \
    -Dcall-caching.enabled=false \
    -Dsystem.job-shell=/bin/sh \
    -jar "${CROMWELL_BUILD_CROMWELL_JAR}" \
    server &

CROMWELL_PID=$!

sleep 30

java \
    -Xmx2g \
    -Dbackend.providers.Local.config.concurrent-job-limit="${CROMWELL_BUILD_CWL_TEST_PARALLELISM}" \
    -jar "${CROMWELL_BUILD_CROMWELL_JAR}" \
    run "${CROMWELL_BUILD_CWL_TEST_WDL}" \
    -i "${CROMWELL_BUILD_CWL_TEST_INPUTS}"

cd "${CROMWELL_BUILD_ROOT_DIRECTORY}"

cromwell::build::generate_code_coverage
