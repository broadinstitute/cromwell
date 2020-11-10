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

# "space" and "scatter" tests are disabled because hey intermittently fail on AWS
# https://broadworkbench.atlassian.net/browse/BA-6152
# NOTE! For some reason, exclusions must all be lower case even if the test name is a mixture of upper and lower
cromwell::build::run_centaur \
    -p 100 \
    -n dummy

cromwell::build::generate_code_coverage
