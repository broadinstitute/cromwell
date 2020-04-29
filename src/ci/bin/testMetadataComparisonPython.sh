#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment
cromwell::build::install_python_scripts_dependencies

build_root="$( dirname "${BASH_SOURCE[0]}" )/../../.."
python3 -m unittest discover -v "${build_root}/scripts/metadata_comparison"

cromwell::build::publish_artifacts
