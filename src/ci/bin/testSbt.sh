#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_unit.inc.sh" || source test_unit.inc.sh

cromwell::build::setup_common_environment

cromwell::build::unit::setup_scale_factor

cromwell::build::unit::setup_exclude_tags

CROMWELL_SBT_TEST_EXCLUDE_TAGS="${CROMWELL_BUILD_UNIT_EXCLUDE_TAGS}" \
CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR="${CROMWELL_BUILD_UNIT_SPAN_SCALE_FACTOR}" \
cromwell:build::run_sbt_test

cromwell::build::generate_code_coverage

cromwell::build::publish_artifacts
