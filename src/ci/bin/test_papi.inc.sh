#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of common Papi functions for use in other scripts.
#
# Functions:
#
#   - cromwell::build::papi::*
#     Functions for use in other Papi scripts
#
#   - cromwell::private::papi::papi::*
#     Functions for use only within this file by cromwell::build::papi::* functions
#

cromwell::build::papi::setup_papi_environment() {
    CROMWELL_BUILD_PAPI_AUTH_JSON="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-service-account.json"
    CROMWELL_BUILD_PAPI_CLIENT_EMAIL="$(jq --exit-status --raw-output .client_email "${CROMWELL_BUILD_PAPI_AUTH_JSON}")"
    CROMWELL_BUILD_PAPI_PROJECT_ID="$(jq --exit-status --raw-output .project_id "${CROMWELL_BUILD_PAPI_AUTH_JSON}")"
    CROMWELL_BUILD_PAPI_GCR_IMAGES="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_papi_gcloud_images_temporary.$$"
    CROMWELL_BUILD_PAPI_CLOUDSDK_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_papi_gcloud_config.$$"

    export CROMWELL_BUILD_PAPI_AUTH_JSON
    export CROMWELL_BUILD_PAPI_CLIENT_EMAIL
    export CROMWELL_BUILD_PAPI_CLIENT_EMAIL_ORIGINAL
    export CROMWELL_BUILD_PAPI_CLOUDSDK_CONFIG
    export CROMWELL_BUILD_PAPI_GCR_IMAGES
    export CROMWELL_BUILD_PAPI_PROJECT_ID

    if [[ "${CROMWELL_BUILD_PROVIDER}" == "${CROMWELL_BUILD_PROVIDER_TRAVIS}" ]]; then
        cromwell::private::papi::install_gcloud
    fi

    # All `gcloud` commands should use this configuration directory.
    # https://stackoverflow.com/questions/34883810/how-to-authenticate-google-apis-with-different-service-account-credentials
    # https://github.com/googleapis/google-auth-library-java/issues/58
    export CLOUDSDK_CONFIG="${CROMWELL_BUILD_PAPI_CLOUDSDK_CONFIG}"

    cromwell::build::add_exit_function cromwell::private::papi::teardown_papi_environment

    gcloud auth activate-service-account --key-file="${CROMWELL_BUILD_PAPI_AUTH_JSON}"
    gcloud config set account "${CROMWELL_BUILD_PAPI_CLIENT_EMAIL}"
    gcloud config set project "${CROMWELL_BUILD_PAPI_PROJECT_ID}"

    if command -v docker; then
        # Upload images built from this commit
        gcloud auth configure-docker --quiet
        CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS="gcr.io/${CROMWELL_BUILD_PAPI_PROJECT_ID}/cromwell-drs-localizer:${CROMWELL_BUILD_CENTAUR_DOCKER_TAG}"
        cromwell::private::papi::gcr_image_push cromwell-drs-localizer "${CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS}"
    else
        # Just use the default images
        CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS="broadinstitute/cromwell-drs-localizer:45-d46ff9f"
    fi

    export CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS
}

cromwell::private::papi::teardown_papi_environment() {
    cromwell::build::delete_docker_images cromwell::private::papi::gcr_image_delete "${CROMWELL_BUILD_PAPI_GCR_IMAGES}"
}

cromwell::private::papi::install_gcloud() {
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
    sudo apt-get install -y apt-transport-https ca-certificates
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
    sudo apt-get update
    sudo apt-get install -y google-cloud-sdk
}

cromwell::private::papi::gcr_image_push() {
    local executable_name
    local docker_image

    executable_name="${1:?gcr_image_push called without an executable_name}"
    docker_image="${2:?gcr_image_push called without an docker_image}"
    shift
    shift

    echo "${docker_image}" >> "${CROMWELL_BUILD_PAPI_GCR_IMAGES}"

    sbt \
        "set \`${executable_name}\`/docker/imageNames := List(ImageName(\"${docker_image}\"))" \
        "${executable_name}/dockerBuildAndPush"
}

cromwell::private::papi::gcr_image_delete() {
    local docker_image_name
    docker_image_name="${1:?gcr_image_delete called without a docker_image_name}"
    shift
    gcloud container images delete "${docker_image_name}" --force-delete-tags --quiet
}
