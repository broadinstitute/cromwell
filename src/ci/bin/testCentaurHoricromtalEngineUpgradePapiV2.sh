#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
export CROMWELL_BUILD_REQUIRES_PULL_REQUEST=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_papi.inc.sh" || source test_papi.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

cromwell::build::papi::setup_papi_environment

prior_version=$(cromwell::private::calculate_prior_version_tag)
export TEST_CROMWELL_PRIOR_VERSION_TAG="${prior_version}"
WOULD_BE_PRIOR_VERSION_CONF="papi_v2_${prior_version}_application.conf"
if [[ -f "$CROMWELL_BUILD_RESOURCES_DIRECTORY/$WOULD_BE_PRIOR_VERSION_CONF" ]]; then
    export TEST_CROMWELL_PRIOR_VERSION_CONF="$WOULD_BE_PRIOR_VERSION_CONF"
else
    export TEST_CROMWELL_PRIOR_VERSION_CONF="papi_v2_application.conf"
fi
# This is the Docker tag that will be applied to the Docker image that is created for the code being built. This image
# will *not* be pushed to Docker Hub or any other repo, it only lives local to the build.
export TEST_CROMWELL_TAG=just-testing-horicromtal
export TEST_CROMWELL_CONF="papi_v2_application.conf"

GOOGLE_AUTH_MODE="service-account"
GOOGLE_REFRESH_TOKEN_PATH="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_refresh_token.txt"

# Export variables used in conf files
export GOOGLE_AUTH_MODE
export GOOGLE_REFRESH_TOKEN_PATH
export TEST_CROMWELL_COMPOSE_FILE="${CROMWELL_BUILD_ROOT_DIRECTORY}/scripts/docker-compose-mysql/docker-compose-horicromtal.yml"

# Copy rendered files
mkdir -p "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"
cp \
    "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2_usa.options.json" \
    "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/private_docker_papi_v2_usa.options" \
    "${TEST_CROMWELL_COMPOSE_FILE}" \
    "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"

docker image ls -q broadinstitute/cromwell:"${TEST_CROMWELL_TAG}" | grep . || \
CROMWELL_SBT_DOCKER_TAGS="${TEST_CROMWELL_TAG}" sbt server/docker

cromwell::build::run_centaur \
    -s "centaur.EngineUpgradeTestCaseSpec" \
    -e localdockertest
