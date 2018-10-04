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

cromwell::build::run_centaur \
    -p 100 \
    -e localdockertest \
    -e gpu_on_papi \
    -d "${CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY}"

cromwell::build::generate_code_coverage
