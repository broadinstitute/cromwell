#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test_gcpbatch.inc.sh" || source test_gcpbatch.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::batch::setup_batch_centaur_environment

cromwell::build::assemble_jars

# Split out of testCentaurGcpBatchRestart.sh, which excludes this test. Restart tests are
# serialized because they kill and restart the shared Cromwell, and this one runs about as
# long as the other three combined.
cromwell::build::run_centaur \
    -p 100 \
    -i call_cache_cha_cha_batch \

cromwell::build::generate_code_coverage
