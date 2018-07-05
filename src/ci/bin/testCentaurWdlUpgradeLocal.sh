#!/usr/bin/env bash

set -e
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

# UpgradeTestCaseSpec takes a selection of ordinary draft-2 test cases (tagged as "upgrade"), runs
# them through the draft-2 to 1.0 upgrade script in Womtool, and runs them against local backend.
centaur/test_cromwell.sh \
    -j "${CROMWELL_BUILD_JAR}" \
    -c "${CROMWELL_BUILD_SCRIPTS_RESOURCES}/local_application.conf" \
    -g \
    -s "centaur.UpgradeTestCaseSpec"
