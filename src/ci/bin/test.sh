#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# Make sure our hooks are configured before building:
git config core.hooksPath | grep -q 'hooks/'

cromwell::build::exec_test_script
