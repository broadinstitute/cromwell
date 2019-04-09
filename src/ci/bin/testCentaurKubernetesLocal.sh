#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# Setting these variables should cause the associated config values to be rendered into centaur_application_horicromtal.conf
# There should probably be more indirections in CI scripts but that can wait.
export TEST_CROMWELL_TAG=just-testing-horicromtal
export TEST_CROMWELL_CONF=horicromtal_application.conf
export CROMWELL_BUILD_MYSQL_USERNAME=travis

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

GOOGLE_AUTH_MODE="service-account"
GOOGLE_REFRESH_TOKEN_PATH="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_refresh_token.txt"

# Export variables used in conf files
export GOOGLE_AUTH_MODE
export GOOGLE_REFRESH_TOKEN_PATH
export TEST_CROMWELL_COMPOSE_FILE="${CROMWELL_BUILD_ROOT_DIRECTORY}/scripts/docker-compose-mysql/docker-compose-horicromtal.yml"

# Copy rendered files
mkdir -p "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"
cp \
    "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/private_docker_papi_v2_usa.options" \
    "${TEST_CROMWELL_COMPOSE_FILE}" \
    "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"

GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-centaur-service-account.json"
gcloud auth activate-service-account --key-file=${GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON}

GOOGLE_KUBERNETES_CLUSTER_NAME=centaur-cluster-${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_NUMBER:-$RANDOM}
GOOGLE_PROJECT=$(docker run --rm -i stedolan/jq:latest < $GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON -r .project_id)

gcloud --project $GOOGLE_PROJECT container clusters create --region us-central $GOOGLE_KUBERNETES_CLUSTER_NAME --num-nodes=3

docker image ls -q broadinstitute/cromwell:"${TEST_CROMWELL_TAG}" | grep . || \
CROMWELL_SBT_DOCKER_TAGS="${TEST_CROMWELL_TAG}" sbt server/docker

cromwell::build::run_centaur \
    -p 100 \
    -e localdockertest \
    "${CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS:-""}" \
    -d "${CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY}"

cromwell::build::generate_code_coverage
