#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_singularity.inc.sh" || source test_singularity.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_slurm.inc.sh" || source test_slurm.inc.sh

cromwell::build::setup_common_environment

cromwell::build::slurm::setup_slurm_environment

cromwell::build::singularity::setup_singularity_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

cromwell::build::run_centaur \
    -e call_cache_capoeira_local \
    -e draft3_call_cache_capoeira_local \
    -e non_root_default_user \
    -e non_root_specified_user \

cromwell::build::generate_code_coverage
