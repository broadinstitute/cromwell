#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of common unit testing functions for use in other scripts.
#
# Functions:
#
#   - cromwell::build::unit::*
#     Functions for use in other unit testing scripts
#
#   - cromwell::private::unit::*
#     Functions for use only within this file by cromwell::build::unit::* functions
#

cromwell::build::unit::setup_scale_factor() {
    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
            CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR=2
            ;;
        "${CROMWELL_BUILD_PROVIDER_JENKINS}")
            CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR=10
            ;;
        *)
            CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR=1
            ;;
    esac
    export CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR
}
