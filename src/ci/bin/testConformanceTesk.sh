#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

if [[ "${CROMWELL_BUILD_EVENT}" != "cron" ]]; then
    echo "TESK integration tests run only as a CRON JOB"
    exit 0
fi

export CROMWELL_BUILD_REQUIRES_SECURE=true

cromwell::build::setup_common_environment
cromwell::build::setup_conformance_environment
cromwell::build::assemble_jars

CENTAUR_CWL_RUNNER_MODE="tesk"
CENTAUR_CWL_JAVA_ARGS="-Dconfig.file=${CROMWELL_BUILD_RESOURCES_DIRECTORY}/ftp_centaur_cwl_runner.conf"
TESK_INPUT_FTP_PREFIX=ftp://ftp.hexdump.org/centaur-cwl-conformance/cwl-inputs/
CROMWELL_BUILD_CWL_TEST_PARALLELISM=8

# Export variables used in conf files and commands
export CENTAUR_CWL_RUNNER_MODE
export CENTAUR_CWL_JAVA_ARGS
export TESK_INPUT_FTP_PREFIX
export CROMWELL_BUILD_CWL_TEST_PARALLELISM

shutdown_cromwell() {
    if [[ -n "${CROMWELL_PID+set}" ]]; then
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
