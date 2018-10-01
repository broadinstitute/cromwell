#!/usr/bin/env bash

set -e
export CROMWELL_BUILD_SUPPORTS_CRON=true
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars


# Installing the AWS CLI
pip install awscli --upgrade --user
export AWS_SHARED_CREDENTIALS_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_credentials
export AWS_CONFIG_FILE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}"/aws_config

# The following tests are skipped:
#
# -i singlesample.aws \
#   Fails as of 2018-09-25 due to an over-eager .stripMargin

cromwell::build::run_centaur \
    -i haplotypcaller.aws

cromwell::build::generate_code_coverage
