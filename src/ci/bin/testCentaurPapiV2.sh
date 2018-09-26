#!/usr/bin/env bash

set -e
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

GOOGLE_AUTH_MODE="service-account"
GOOGLE_REFRESH_TOKEN_PATH="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_refresh_token.txt"

# Export variables used in conf files
export GOOGLE_AUTH_MODE
export GOOGLE_REFRESH_TOKEN_PATH

# Copy rendered files
mkdir -p "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"
cp \
    "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/private_docker_papi_v2_usa.options" \
    "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"

# Excluded tests:
# docker_hash_dockerhub_private: https://github.com/broadinstitute/cromwell/issues/3587

cromwell::build::run_centaur \
    -p 100 \
    -e localdockertest \
    -d "${CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY}"

cromwell::build::generate_code_coverage
