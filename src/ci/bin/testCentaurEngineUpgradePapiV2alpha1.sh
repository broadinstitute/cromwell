#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
export CROMWELL_BUILD_REQUIRES_PRIOR_VERSION=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_papi.inc.sh" || source test_papi.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::papi::setup_papi_centaur_environment

cromwell::build::assemble_jars

cromwell::build::run_centaur \
    -s "centaur.EngineUpgradeTestCaseSpec" \
    -e localdockertest \

cromwell::build::generate_code_coverage
