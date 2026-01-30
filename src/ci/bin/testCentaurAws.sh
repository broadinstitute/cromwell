#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

export AWS_SHARED_CREDENTIALS_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_credentials
export AWS_CONFIG_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_config

# Get AWS credentials from Vault
export AWS_ACCESS_KEY=$(vault read -field=access_key secret/dsde/cromwell/common/cromwell-aws)
export AWS_SECRET_KEY=$(vault read -field=secret_key secret/dsde/cromwell/common/cromwell-aws)

# TODO turn most tests back on once we resolve timeouts
# TODO (AN-710) Add back some of these tests (space, scatter, docker_hash_dockerhub, awswdlresultscopying etc.)
# TODO (AN-710) tests that depend on continueOnReturnCode tests are failing:
# (exit, valid_return_codes_and_continue_on_return_code, return_codes, globbingBehavior, failures.terminal_status)
# TODO (AN-794) support GPU tests in AWS job queue (enables test gpu_required_and_requested)
cromwell::build::run_centaur \
    -p 100 \
    -e localdockertest \
    -e abort.scheduled_abort \
    -e relative_output_paths \
    -e relative_output_paths_colliding \
    -e standard_output_paths_colliding_prevented \
    -e restart \
    -e space \
    -e scatter \
    -e runtwiceexpectingcallcaching \
    -e papi_v2alpha1_gcsa \
    -e docker_hash_dockerhub \
    -e awswdlresultscopying \
    -e awswdlresultscopyingrelative \
    -e exit \
    -e default_runtime_attributes \
    -e valid_return_codes_and_continue_on_return_code \
    -e dont_cache_to_failed_jobs \
    -e continue_on_return_code \
    -e return_codes \
    -e globbingbehavior \
    -e cachewithinwf \
    -e failures.terminal_status \
    -e bad_file_string \
    -e awsbatch_labels \
    -e gpu_required_and_requested


cromwell::build::generate_code_coverage

cromwell::build::print_workflow_statistics
