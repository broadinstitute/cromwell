#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

if [[ "${CROMWELL_BUILD_EVENT}" != "cron" ]]; then
    echo "TESK integration tests run only as a CRON JOB"
    exit 0
fi

cromwell::build::setup_common_environment

cromwell::build::setup_conformance_environment

cromwell::build::assemble_jars

CENTAUR_CWL_JAVA_ARGS="-Dconfig.file=${CROMWELL_BUILD_RESOURCES_DIRECTORY}/ftp_centaur_cwl_runner.conf"
TESK_INPUT_FTP_PREFIX=ftp://ftp.hexdump.org/centaur-cwl-conformance/cwl-inputs/
CROMWELL_BUILD_CWL_TEST_PARALLELISM=8

export CENTAUR_CWL_JAVA_ARGS
export TESK_INPUT_FTP_PREFIX
export CROMWELL_BUILD_CWL_TEST_PARALLELISM

cromwell::build::run_conformance

cromwell::build::generate_code_coverage
