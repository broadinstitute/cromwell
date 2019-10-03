#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

# WdlUpgradeTestCaseSpec takes a selection of ordinary draft-2 test cases (tagged as "upgrade"), runs
# them through the draft-2 to 1.0 upgrade script in Womtool, and runs them against local backend.
cromwell::build::run_centaur \
    -s "centaur.WdlUpgradeTestCaseSpec"

cromwell::build::generate_code_coverage
