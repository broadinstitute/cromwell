#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
set -x
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::exec_test_script
