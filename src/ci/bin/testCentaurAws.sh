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
    -e smartseq2singlesample \
    -e arrays \
    -e haplotypecaller \
    -e jointdiscovery \
    -e mutect2 \
    -e singlesample \
    -e singlesample_production \
    -e cnv_somatic_pair \
    -e cnv_somatic_panel \
    -e localdockertest \
    -e inline_file \
    -e iwdr_input_string_function \
    -e non_root_default_user \
    -e cwl_interpolated_strings \
    -e non_root_specified_user \
    -e space \
    -e draft3_optional_input_from_scatter \
    -e iwdr_input_string \
    -e cwl_cache_between_workflows \
    -e abort.scheduled_abort \
    -e cwl_cache_within_workflow \
    -e inline_file_custom_entryname \
    -e relative_output_paths \
    -e relative_output_paths_colliding \
    -e standard_output_paths_colliding_prevented \

cromwell::build::generate_code_coverage
