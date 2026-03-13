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

cromwell::build::run_centaur \
    -p 100 \
    -i awsbatch_labels \

cromwell::build::generate_code_coverage

cromwell::build::print_workflow_statistics
