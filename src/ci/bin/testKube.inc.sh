#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of functions for use in the Kubernetes test script.

# Takes a single string argument and `echo`s a possibly modified version of that argument with non-alphanumeric
# characters converted to dashes. TODO: restrict the initial character as necessary
cromwell::kube::google_safe_name() {
  echo -n "$1" | tr -c '[[:digit:][:alpha:]]' '-'
}

# Creates a Google friendly identifier name specific to this build based on its single argument.
cromwell::kube::centaur_gke_name() {
  local prefix="centaur-gke"
  local build_name="$(cromwell::kube::google_safe_name ${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_NUMBER:-$RANDOM})"
  local arg=$1
  echo -n "${prefix}-${arg}-${build_name}" | tr -c '[[:digit:][:alpha:]]' '-'
}

# Run a specified command after activating the specified service account.
#
# Usage: cromwell::kube::gcloud_run_as_service_account command service_account.json
cromwell::kube::gcloud_run_as_service_account() {
  local command="$1"
  local service_json="$2"
  local DOCKER_ETC_PATH=/usr/share/etc
  docker run -v "$CROMWELL_BUILD_RESOURCES_DIRECTORY:$DOCKER_ETC_PATH" -e DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/${service_json} && $command "
}
