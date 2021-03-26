#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_unit.inc.sh" || source test_unit.inc.sh

cromwell::build::setup_common_environment

cromwell::build::unit::setup_unit_environment

# Give unit tests more memory via SBT_OPTS leaving .sbtopts with the default for all other integration tests.
# As there should be no MySQL, PostgreSQL, or other DBMS containers running in this CI we should have more memory
# available. Hopefully this extra memory will help stabalize non-dbms unit tests that sometimes stress memory and CPU.
# For more information on testing and memory see also: https://olegych.github.io/blog/sbt-fork.html

CROMWELL_SBT_TEST_EXCLUDE_TAGS="${CROMWELL_BUILD_UNIT_EXCLUDE_TAGS}" \
CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR="${CROMWELL_BUILD_UNIT_SPAN_SCALE_FACTOR}" \
SBT_OPTS="${SBT_OPTS:-} ${CROMWELL_BUILD_UNIT_SBT_MEMORY}" \
cromwell:build::run_sbt_test

cromwell::build::generate_code_coverage

cromwell::build::publish_artifacts
