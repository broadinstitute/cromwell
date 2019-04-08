#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true

# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_kube.inc.sh" || source test_kube.inc.sh

# Comment image build out during change lock debugging
#DUMMY_IMAGE="just-testing"
#KUBE_CROMWELL_IMAGE=$(cromwell::kube::generate_gcr_tag)
#CROMWELL_SBT_DOCKER_TAGS=${DUMMY_IMAGE} sbt server/docker
#cromwell::kube::gcr_login
#cromwell::kube::push_to_gcr ${DUMMY_IMAGE} ${KUBE_CROMWELL_IMAGE}
KUBE_CROMWELL_IMAGE=broadinstitute/cromwell:41


KUBE_CLOUDSQL_INSTANCE_NAME="$(cromwell::kube::generate_cloud_sql_instance_name)"
cromwell::kube::create_cloud_sql_instance ${KUBE_CLOUDSQL_INSTANCE_NAME}
# Get the connectionName for this newly created instance. This is what the Cloud SQL proxies will need for their -instances parameter.
KUBE_CLOUDSQL_CONNECTION_NAME="$(cromwell::kube::connection_name_for_cloud_sql_instance ${KUBE_CLOUDSQL_INSTANCE_NAME})"
echo "Cloud SQL connection name is $KUBE_CLOUDSQL_CONNECTION_NAME"


KUBE_CLUSTER_NAME=$(cromwell::kube::generate_gke_cluster_name)
cromwell::kube::create_gke_cluster ${KUBE_CLUSTER_NAME}
cromwell::kube::create_deployment_files
cromwell::kube::create_secrets

cromwell::kube::start_cromwell

cromwell::kube::get_load_balancer_ip "ip.txt"
KUBE_LOAD_BALANCER_IP=$(cat "ip.txt")

# Setting these variables should cause the associated config values to be rendered into centaur_application_horicromtal_docker_compose.conf
# There should probably be more indirections in CI scripts but that can wait.
export TEST_CROMWELL_CONF=horicromtal_application.conf
export CROMWELL_BUILD_MYSQL_USERNAME=travis


#cromwell::build::assemble_jars

GOOGLE_AUTH_MODE="service-account"
GOOGLE_REFRESH_TOKEN_PATH="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_refresh_token.txt"

# Export variables used in conf files
export GOOGLE_AUTH_MODE
export GOOGLE_REFRESH_TOKEN_PATH
export TEST_CROMWELL_COMPOSE_FILE="${CROMWELL_BUILD_ROOT_DIRECTORY}/scripts/docker-compose-mysql/docker-compose-horicromtal.yml"

# this is probably vestigial and should just be deleted but not 100% sure
# Copy rendered files
#mkdir -p "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"
#cp \
#    "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/private_docker_papi_v2_usa.options" \
#    "${TEST_CROMWELL_COMPOSE_FILE}" \
#    "${CROMWELL_BUILD_CENTAUR_TEST_RENDERED}"

# TODO Move this to the "cleanup" section of the script once there is also a "do real work" section.
# cromwell:kube::destroy_gke_cluster ${KUBE_CLUSTER_NAME}

## TODO Move this to the "cleanup" section of the script once there is also a "do real work" section.
#cromwell::kube::destroy_cloud_sql_instance ${KUBE_CLOUDSQL_INSTANCE_NAME}

# Temporarily turning GCR push off as it's unnecessary for testing Cloud SQL proxy stuff.
#GCR_TAG=$(cromwell::kube::generate_gcr_tag)
#cromwell::kube::tag_for_gcr ${DOCKER_IMAGE} ${GCR_TAG}
#cromwell::kube::gcr_login
#cromwell::kube::push_to_gcr ${GCR_TAG}
# TODO Move this to the "cleanup" section of the script once there is also a "do real work" section.
#cromwell::kube::delete_from_gcr ${GCR_TAG}

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

