#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

if [ "${CROMWELL_BUILD_PROVIDER}" = "${CROMWELL_BUILD_PROVIDER_TRAVIS}" ] && [ -n "${TRAVIS_PULL_REQUEST_BRANCH}" ]; then

  cromwell::build::setup_common_environment

  cromwell::build::setup_centaur_environment

  cromwell::build::assemble_jars

  cromwell::build::run_centaur \
      -p 100 \
      -e localdockertest \
      -e gpu_on_papi \

fi
