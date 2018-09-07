#!/usr/bin/env bash

set -e
export CROMWELL_BUILD_SUPPORTS_CRON=true
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

GOOGLE_AUTH_MODE="service-account"
GOOGLE_REFRESH_TOKEN_PATH="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_refresh_token.txt"
GOOGLE_SERVICE_ACCOUNT_JSON="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-service-account.json"

# Export variables used in conf files
export GOOGLE_AUTH_MODE
export GOOGLE_REFRESH_TOKEN_PATH
export GOOGLE_SERVICE_ACCOUNT_JSON

# pass integration directory to the inputs json otherwise remove it from the inputs file
INTEGRATION_TESTS=()
if [ "${CROMWELL_BUILD_IS_CRON}" = "true" ]; then
    INTEGRATION_TESTS=(-d "${CROMWELL_BUILD_CENTAUR_INTEGRATION_TESTS}")
    # Increase concurrent job limit to get tests to finish under three hours.
    # Increase read_lines limit because of read_lines call on hg38.even.handcurated.20k.intervals.
    CENTAUR_READ_LINES_LIMIT=512000
    export CENTAUR_READ_LINES_LIMIT
fi

centaur/test_cromwell.sh \
    -j "${CROMWELL_BUILD_JAR}" \
    -c "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v1_application.conf" \
    -n "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application.conf" \
    -p 100 \
    -g \
    -e localdockertest \
    -e gpu_on_papi \
    "${INTEGRATION_TESTS[@]}"

cromwell::build::generate_code_coverage
