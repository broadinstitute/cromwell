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

cromwell::private::papi::setup_papi_gcloud() {
    CROMWELL_BUILD_PAPI_AUTH_JSON="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-service-account.json"
    CROMWELL_BUILD_PAPI_CLIENT_EMAIL="$(jq --exit-status --raw-output .client_email "${CROMWELL_BUILD_PAPI_AUTH_JSON}")"
    CROMWELL_BUILD_PAPI_PROJECT_ID="$(jq --exit-status --raw-output .project_id "${CROMWELL_BUILD_PAPI_AUTH_JSON}")"
    CROMWELL_BUILD_PAPI_GCR_IMAGES="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_papi_gcloud_images_temporary.$$"
    CROMWELL_BUILD_PAPI_CLOUDSDK_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_papi_gcloud_config.$$"

    export CROMWELL_BUILD_PAPI_AUTH_JSON
    export CROMWELL_BUILD_PAPI_CLIENT_EMAIL
    export CROMWELL_BUILD_PAPI_CLOUDSDK_CONFIG
    export CROMWELL_BUILD_PAPI_GCR_IMAGES
    export CROMWELL_BUILD_PAPI_PROJECT_ID

    # All `gcloud` commands should use this configuration directory.
    # https://stackoverflow.com/questions/34883810/how-to-authenticate-google-apis-with-different-service-account-credentials
    # https://github.com/googleapis/google-auth-library-java/issues/58
    export CLOUDSDK_CONFIG="${CROMWELL_BUILD_PAPI_CLOUDSDK_CONFIG}"

    cromwell::build::add_exit_function cromwell::private::papi::teardown_papi_gcloud

    gcloud auth activate-service-account --key-file="${CROMWELL_BUILD_PAPI_AUTH_JSON}"
    gcloud config set account "${CROMWELL_BUILD_PAPI_CLIENT_EMAIL}"
    gcloud config set project "${CROMWELL_BUILD_PAPI_PROJECT_ID}"
}

cromwell::private::papi::teardown_papi_gcloud() {
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

    cromwell::build::build_docker_image "${executable_name}" "${docker_image}"
    echo "${docker_image}" >> "${CROMWELL_BUILD_PAPI_GCR_IMAGES}"
    # Use cat to quiet docker: https://github.com/moby/moby/issues/36655#issuecomment-375136087
    docker push "${docker_image}" | cat
}

cromwell::private::papi::gcr_image_delete() {
    local docker_image_name
    docker_image_name="${1:?gcr_image_delete called without a docker_image_name}"
    shift
    gcloud container images delete "${docker_image_name}" --force-delete-tags --quiet
}

cromwell::private::papi::setup_papi_gcr() {
    # Build a DOS/DRS localizer image from source, or for faster local debugging use an already provided image
    if [[ -n "${CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS:+set}" ]]; then
        # If CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS is already set then use that image
        echo "Using CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS='${CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS}'"
    elif command -v docker; then
        # Upload images built from this commit
        gcloud auth configure-docker --quiet
        CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS="gcr.io/${CROMWELL_BUILD_PAPI_PROJECT_ID}/cromwell-drs-localizer:${CROMWELL_BUILD_DOCKER_TAG}-papi"
        cromwell::private::papi::gcr_image_push cromwell-drs-localizer "${CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS}"
        export CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS
    else
        echo "Error: BA-6546 The environment variable CROMWELL_BUILD_PAPI_DOCKER_IMAGE_DRS must be set/export pointing to a valid docker image" >&2
        exit 1
    fi
}

cromwell::private::papi::setup_papi_service_account() {
    CROMWELL_BUILD_PAPI_AUTH_MODE="service-account"
    export CROMWELL_BUILD_PAPI_AUTH_MODE
}

cromwell::private::papi::setup_papi_endpoint_url() {
    if [[ "${CROMWELL_BUILD_TYPE}" == "centaurPapiV2" ]]; then
        CROMWELL_BUILD_PAPI_ENDPOINT_URL="https://lifesciences.googleapis.com/"
    else
        CROMWELL_BUILD_PAPI_ENDPOINT_URL="https://genomics.googleapis.com/"
    fi

    export CROMWELL_BUILD_PAPI_ENDPOINT_URL
}

cromwell::build::papi::setup_papi_centaur_environment() {
    cromwell::private::papi::setup_papi_gcloud
    cromwell::private::papi::setup_papi_gcr
    cromwell::private::papi::setup_papi_service_account
    cromwell::private::papi::setup_papi_endpoint_url
}
