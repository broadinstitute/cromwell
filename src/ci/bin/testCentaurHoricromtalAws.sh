#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_aws.inc.sh" || source test_aws.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

#cromwell::build::batch::setup_aws_centaur_environment

cromwell::build::assemble_jars

cromwell::build::build_cromwell_docker

export AWS_SHARED_CREDENTIALS_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_credentials
export AWS_CONFIG_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_config


cromwell::build::run_centaur \
    -p 500 \
    -e localdockertest \
    -e relative_output_paths \
    -e relative_output_paths_colliding \
    -e standard_output_paths_colliding_prevented \
    -e papi_v2alpha1_gcsa \
    -e restart \
    -e lots_of_inputs_papiv2 \
    # Due to Centaur changes in WX-1629, lots_of_inputs_papiv2 times out due to the large number of inputs. Instead,
    # we run lots_of_inputs, which tests 400 inputs instead of the 10,000 in lots_of_inputs_papiv2

cromwell::build::generate_code_coverage

cromwell::build::print_workflow_statistics
