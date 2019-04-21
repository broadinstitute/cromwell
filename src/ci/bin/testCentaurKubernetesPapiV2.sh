#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/testKube.inc.sh" || source testKube.inc.sh

# Setting these variables should cause the associated config values to be rendered into centaur_application_horicromtal.conf
# There should probably be more indirections in CI scripts but that can wait.
export TEST_CROMWELL_TAG=just-testing-horicromtal
export TEST_CROMWELL_CONF=horicromtal_application.conf
export CROMWELL_BUILD_MYSQL_USERNAME=travis


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

# Temporarily turning GKE cluster stuff off as it's unnecessary for testing Cloud SQL stuff.
#
# KUBE_CLUSTER_NAME=$(cromwell::kube::generate_gke_cluster_name)
# cromwell::kube::create_gke_cluster ${KUBE_CLUSTER_NAME}
# TODO Move this to the "cleanup" section of the script once there is also a "do real work" section.
# cromwell:kube::destroy_gke_cluster ${KUBE_CLUSTER_NAME}

# Temporarily turning off Cloud SQL stuff as it's unnecessary for testing GCR stuff.
#KUBE_CLOUDSQL_INSTANCE_NAME="$(cromwell::kube::generate_cloud_sql_instance_name)"
#echo "Cloud SQL instance name is $KUBE_CLOUDSQL_INSTANCE_NAME"
#
#cromwell::kube::create_cloud_sql_instance ${KUBE_CLOUDSQL_INSTANCE_NAME}
## Get the connectionName for this newly created instance. This is what the Cloud SQL proxies will need for their -instances parameter.
#KUBE_CLOUDSQL_CONNECTION_NAME="$(cromwell::kube::connection_name_for_cloud_sql_instance ${KUBE_CLOUDSQL_INSTANCE_NAME})"
#echo "Cloud SQL connection name is $KUBE_CLOUDSQL_CONNECTION_NAME"
#
## TODO Move this to the "cleanup" section of the script once there is also a "do real work" section.
#cromwell::kube::destroy_cloud_sql_instance ${KUBE_CLOUDSQL_INSTANCE_NAME}

# Just use Cromwell 40 for testing pushes to GCR. IRL this script would obvs build the image and push that.
docker pull broadinstitute/cromwell:40
KUBE_GCR_TAG=$(cromwell::kube::generate_gcr_tag)
docker tag broadinstitute/cromwell:40 ${KUBE_GCR_IMAGE_TAG}
docker push ${KUBE_GCR_IMAGE_TAG}
# TODO Move this to the "cleanup" section of the script once there is also a "do real work" section.
cromwell::kube::delete_gcr_image ${KUBE_GCR_IMAGE_TAG}

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
