#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_gcpbatch.inc.sh" || source test_gcpbatch.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::batch::setup_batch_centaur_environment

cromwell::build::assemble_jars

cromwell::build::run_centaur \
    -p 100 \
    -i restart \

cromwell::build::generate_code_coverage

cromwell::build::print_workflow_statistics
