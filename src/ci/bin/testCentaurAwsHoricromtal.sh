#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

cromwell::build::build_cromwell_docker

export AWS_SHARED_CREDENTIALS_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_credentials
export AWS_CONFIG_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_config

# Get AWS credentials from Vault
export AWS_ACCESS_KEY=$(vault read -field=access_key secret/dsde/cromwell/common/cromwell-aws)
export AWS_SECRET_KEY=$(vault read -field=secret_key secret/dsde/cromwell/common/cromwell-aws)

# TODO turn most tests back on once we resolve timeouts
# TODO (AN-710) Add back some of these tests (space, scatter, docker_hash_dockerhub, awswdlresultscopying etc.)
cromwell::build::run_centaur \
     -i hello \
     -i mutect2.aws
#    -p 500 \
#    -e localdockertest \
#    -e abort.scheduled_abort \
#    -e relative_output_paths \
#    -e relative_output_paths_colliding \
#    -e standard_output_paths_colliding_prevented \
#    -e restart \
#    -e space \
#    -e scatter \
#    -e runtwiceexpectingcallcaching \
#    -e papi_v2alpha1_gcsa \
#    -e docker_hash_dockerhub \
#    -e awswdlresultscopying \
#    -e awswdlresultscopyingrelative

cromwell::build::generate_code_coverage

cromwell::build::print_workflow_statistics
