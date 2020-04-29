#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_unit.inc.sh" || source test_unit.inc.sh

cromwell::build::setup_common_environment
cromwell::build::install_python_scripts_dependencies

python3 -m unittest discover -v ${BASH_SOURCE%/*}/metadata_comparison

cromwell::build::publish_artifacts
