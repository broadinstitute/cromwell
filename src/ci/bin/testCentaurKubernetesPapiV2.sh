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

#cromwell::build::assemble_jars

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
GOOGLE_ZONE=us-central1-c

KUBE_CLUSTER_NAME=$(centaur_gke_name "cluster")
KUBE_SQL_INSTANCE_NAME=$(centaur_gke_name "cloudsql")
KUBE_CLOUDSQL_PASSWORD="$(docker run --rm -e VAULT_TOKEN=${VAULT_TOKEN} broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/centaur-gke | jq -r '.data.db_pass' | tr -d '\n')"

GOOGLE_PROJECT=$(docker run --rm -i stedolan/jq:latest < $GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON -r .project_id)

# TEMP TURNING THIS OFF TO TEST CLOUDSQL STUFF
# gcloud --project $GOOGLE_PROJECT container clusters create --zone $GOOGLE_ZONE $KUBE_CLUSTER_NAME --num-nodes=3
#WARNING: Accessing a Container Engine cluster requires the kubernetes commandline
#client [kubectl]. To install, run
#  $ gcloud components install kubectl

# Phase 1. Even this is PAPI since Cromwell will be running in a Docker container and trying to run Docker in Docker
#          currently no es bueno.
# - spin up a Cloud SQL. Obtain its coordinates to be able to access it from a Cloud SQL proxy.
#   (I think this might have been why I didn't do Cloud SQL before but who cares I think it's worth it).
gcloud --project $GOOGLE_PROJECT sql instances create --zone $GOOGLE_ZONE --storage-size=10GB $KUBE_SQL_INSTANCE_NAME

# Create a user
gcloud --project $GOOGLE_PROJECT sql users create cromwell --instance $KUBE_SQL_INSTANCE_NAME --password="${KUBE_CLOUDSQL_PASSWORD}"

# - spin up a CloudIP service fronting said MySQL container
# - spin up a uni-Cromwell that talks to said MySQL
# - spin up a LoadBalancer service that fronts this Cromwell
#
# Run the Centaur test suite against this Cromwell service.

# Phase 1. Even this is PAPI since Cromwell will be running in a Docker container and trying to run Docker in Docker
#          currently no es bueno.
# - spin up a MySQL container expressing a PersistentVolumeClaim
# - spin up a CloudIP service fronting said MySQL container
# - spin up a uni-Cromwell that talks to said MySQL
# - spin up a LoadBalancer service that fronts this Cromwell
#
# Run the Centaur test suite against this Cromwell service.

# Phase 2 same as Phase 1 except separate Cromwells for summarizer, frontend, backend.

# TEMP TURNING THIS OFF TO TEST CLOUDSQL STUFF
# gcloud --project $GOOGLE_PROJECT --quiet container clusters delete $GOOGLE_KUBERNETES_CLUSTER_NAME --zone $GOOGLE_ZONE

#docker image ls -q broadinstitute/cromwell:"${TEST_CROMWELL_TAG}" | grep . || \
#CROMWELL_SBT_DOCKER_TAGS="${TEST_CROMWELL_TAG}" sbt server/docker
#
#cromwell::build::run_centaur \
#    -p 100 \
#    -e localdockertest \
#    "${CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS:-""}" \
#    -d "${CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY}"
#
#cromwell::build::generate_code_coverage

# Creates a build instance specific, Google friendly identifier name based on its sole argument.
centaur_gke_name() {
  local prefix="centaur-gke"
  local build_name="$(cromwell::build::google_safe_name ${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_NUMBER:-$RANDOM})"
  local arg=$1
  echo -n "${prefix}-${arg}-${build_name}" | tr -c '[[:digit:][:alpha:]]' '-'
}
