#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

# After fighting with substitutions vs. environment variables using this wrapper script ONLY sets environment variables

substitutions="_REPO_NAME=${REPO_NAME:-}"
substitutions="${substitutions},_BRANCH_NAME=${BRANCH_NAME:-}"
substitutions="${substitutions},_TAG_NAME=${TAG_NAME:-}"
substitutions="${substitutions},_PR_NUMBER=${PR_NUMBER:-}"
substitutions="${substitutions},_BUILD_NAME=${BUILD_NAME:-}"
substitutions="${substitutions},_BUILD_TYPE=${BUILD_TYPE:-}"
substitutions="${substitutions},_BUILD_HSQLDB=${BUILD_HSQLDB:-}"
substitutions="${substitutions},_BUILD_MARIADB=${BUILD_MARIADB:-}"
substitutions="${substitutions},_BUILD_MYSQL=${BUILD_MYSQL:-}"
substitutions="${substitutions},_BUILD_POSTGRESQL=${BUILD_POSTGRESQL:-}"
substitutions="${substitutions},_BUILD_SQLITE=${BUILD_SQLITE:-}"

# Allow the .git directory to be uploaded by not ignoring it
# https://github.com/GoogleCloudPlatform/cloud-builders/issues/236
gcloud config set gcloudignore/enabled false
# NOTE: This assumes this script is run from the root
gcloud builds submit --config=cloudbuild_single.yaml --substitutions="${substitutions}" .
