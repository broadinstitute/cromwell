#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_papi.inc.sh" || source test_papi.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::papi::setup_papi_centaur_environment

cromwell::build::assemble_jars

cromwell::build::build_cromwell_docker

cromwell::build::run_centaur \
    -p 100 \
    -e localdockertest \
    -e relative_output_paths \
    -e relative_output_paths_colliding \
    -e standard_output_paths_colliding_prevented \

cromwell::build::generate_code_coverage
