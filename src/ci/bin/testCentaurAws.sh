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
    -p 100 \
    -e abort.scheduled_abort \
    -e arrays \
    -e bad_file_string \
    -e bad_output_task \
    -e cnv_somatic_pair \
    -e cnv_somatic_panel \
    -e continue_on_return_code \
    -e cwl_cache_between_workflows \
    -e cwl_cache_within_workflow \
    -e cwl_interpolated_strings \
    -e default_runtime_attributes \
    -e draft3_glob_access \
    -e draft3_globs \
    -e draft3_optional_input_from_scatter \
    -e exit \
    -e failures.terminal_status \
    -e file_name_too_long \
    -e globbingbehavior \
    -e globbingindex \
    -e globbingscatter \
    -e haplotypecaller \
    -e inline_file \
    -e inline_file_custom_entryname \
    -e iwdr_input_string \
    -e iwdr_input_string_function \
    -e jointdiscovery \
    -e localdockertest \
    -e mutect2 \
    -e non_root_default_user \
    -e non_root_specified_user \
    -e relative_output_paths \
    -e relative_output_paths_colliding \
    -e singlesample \
    -e singlesample_production \
    -e smartseq2singlesample \
    -e space \
    -e standard_output_paths_colliding_prevented \

cromwell::build::generate_code_coverage
