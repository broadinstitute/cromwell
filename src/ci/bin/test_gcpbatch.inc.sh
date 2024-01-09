#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of common Gcp Batch functions for use in other scripts.
#
# Functions:
#
#   - cromwell::build::batch::*
#     Functions for use in other Papi scripts
#
#   - cromwell::private::batch::batch::*
#     Functions for use only within this file by cromwell::build::batch::* functions
#

cromwell::private::batch::setup_batch_gcloud() {
    CROMWELL_BUILD_BATCH_AUTH_JSON="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-service-account.json"
    CROMWELL_BUILD_BATCH_CLIENT_EMAIL="$(jq --exit-status --raw-output .client_email "${CROMWELL_BUILD_BATCH_AUTH_JSON}")"
    CROMWELL_BUILD_BATCH_PROJECT_ID="$(jq --exit-status --raw-output .project_id "${CROMWELL_BUILD_BATCH_AUTH_JSON}")"
    CROMWELL_BUILD_BATCH_GCR_IMAGES="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_batch_gcloud_images_temporary.$$"
    CROMWELL_BUILD_BATCH_CLOUDSDK_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_batch_gcloud_config.$$"

    export CROMWELL_BUILD_BATCH_AUTH_JSON
    export CROMWELL_BUILD_BATCH_CLIENT_EMAIL
    export CROMWELL_BUILD_BATCH_CLOUDSDK_CONFIG
    export CROMWELL_BUILD_BATCH_GCR_IMAGES
    export CROMWELL_BUILD_BATCH_PROJECT_ID

    # All `gcloud` commands should use this configuration directory.
    # https://stackoverflow.com/questions/34883810/how-to-authenticate-google-apis-with-different-service-account-credentials
    # https://github.com/googleapis/google-auth-library-java/issues/58
    export CLOUDSDK_CONFIG="${CROMWELL_BUILD_BATCH_CLOUDSDK_CONFIG}"

    cromwell::build::add_exit_function cromwell::private::batch::teardown_batch_gcloud

    gcloud auth activate-service-account --key-file="${CROMWELL_BUILD_BATCH_AUTH_JSON}"
    export GOOGLE_APPLICATION_CREDENTIALS="${CROMWELL_BUILD_BATCH_AUTH_JSON}"
    gcloud config set account "${CROMWELL_BUILD_BATCH_CLIENT_EMAIL}"
    gcloud config set project "${CROMWELL_BUILD_BATCH_PROJECT_ID}"
}

cromwell::private::batch::teardown_batch_gcloud() {
    cromwell::build::delete_docker_images cromwell::private::batch::gcr_image_delete "${CROMWELL_BUILD_BATCH_GCR_IMAGES}"
}

cromwell::private::batch::install_gcloud() {
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
    sudo apt-get install -y apt-transport-https ca-certificates
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
    sudo apt-get update
    sudo apt-get install -y google-cloud-sdk
}

cromwell::private::batch::gcr_image_push() {
    local executable_name
    local docker_image

    executable_name="${1:?gcr_image_push called without an executable_name}"
    docker_image="${2:?gcr_image_push called without an docker_image}"
    shift
    shift

    cromwell::build::build_docker_image "${executable_name}" "${docker_image}"
    echo "${docker_image}" >> "${CROMWELL_BUILD_BATCH_GCR_IMAGES}"
    # Use cat to quiet docker: https://github.com/moby/moby/issues/36655#issuecomment-375136087
    docker push "${docker_image}" | cat
}

cromwell::private::batch::gcr_image_delete() {
    local docker_image_name
    docker_image_name="${1:?gcr_image_delete called without a docker_image_name}"
    shift
    gcloud container images delete "${docker_image_name}" --force-delete-tags --quiet
}

cromwell::private::batch::setup_batch_gcr() {
    # Build a DOS/DRS localizer image from source, or for faster local debugging use an already provided image
    if [[ -n "${CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS:+set}" ]]; then
        # If CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS is already set then use that image
        echo "Using CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS='${CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS}'"
    elif command -v docker; then
        # Upload images built from this commit
        gcloud auth configure-docker --quiet
        CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS="gcr.io/${CROMWELL_BUILD_BATCH_PROJECT_ID}/cromwell-drs-localizer:${CROMWELL_BUILD_DOCKER_TAG}-batch"
        cromwell::private::batch::gcr_image_push cromwell-drs-localizer "${CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS}"
        export CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS
    else
        echo "Error: BA-6546 The environment variable CROMWELL_BUILD_BATCH_DOCKER_IMAGE_DRS must be set/export pointing to a valid docker image" >&2
        exit 1
    fi
}

cromwell::private::batch::setup_batch_service_account() {
    CROMWELL_BUILD_BATCH_AUTH_MODE="service-account"
    export CROMWELL_BUILD_BATCH_AUTH_MODE
}

cromwell::build::batch::setup_batch_centaur_environment() {
    cromwell::private::batch::setup_batch_gcloud
    cromwell::private::batch::setup_batch_service_account
}
