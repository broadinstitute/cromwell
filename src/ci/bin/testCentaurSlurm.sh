#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_slurm.inc.sh" || source test_slurm.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

cromwell::build::slurm::slurm_install

cromwell::build::run_centaur \
    -e call_cache_capoeira_local \
    -e draft3_call_cache_capoeira_local \

cromwell::build::generate_code_coverage
