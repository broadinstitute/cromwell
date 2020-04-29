#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_unit.inc.sh" || source test_unit.inc.sh

echo 1
cromwell::build::setup_common_environment
echo 2
cromwell::build::install_python_scripts_dependencies
echo 3

python3 -m unittest discover -v ${BASH_SOURCE%/*}/metadata_comparison
echo 4

cromwell::build::publish_artifacts
echo 5
