#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_SUPPORTS_CRON=true
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars


# Installing the AWS CLI
cromwell::build::pip_install awscli --upgrade
export AWS_SHARED_CREDENTIALS_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_credentials
export AWS_CONFIG_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_config


cromwell::build::run_centaur \
    -e localdockertest \
    -e inline_file \
    -e exit \
    -e iwdr_input_string_function \
    -e globbingBehavior \
    -e non_root_default_user \
    -e draft3_glob_access \
    -e bad_output_task \
    -e cacheWithinWF \
    -e floating_tags \
    -e cwl_interpolated_strings \
    -e cacheBetweenWF \
    -e failures.terminal_status \
    -e call_cache_hit_prefixes_no_hint \
    -e bad_file_string \
    -e default_runtime_attributes \
    -e non_root_specified_user \
    -e space \
    -e draft3_optional_input_from_scatter \
    -e iwdr_input_string \
    -e globbingindex \
    -e cwl_cache_within_workflow \
    -e continue_on_return_code \
    -e globbingscatter \
    -e inline_file_custom_entryname \
    -e draft3_globs \    
    -e cwl_cache_between_workflows \
    -e abort.scheduled_abort


cromwell::build::generate_code_coverage
