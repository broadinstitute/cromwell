#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

#NB: This function ensures that the correct .conf file is being used by cromwell (among other things).
#Blob storage requires a configuration file tailored for Azure.
cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

cromwell::build::run_centaur

cromwell::build::generate_code_coverage
