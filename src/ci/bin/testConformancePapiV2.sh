#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_conformance_environment

cromwell::build::assemble_jars

GOOGLE_AUTH_MODE="service-account"
PAPI_INPUT_GCS_PREFIX=gs://centaur-cwl-conformance-1f501e3/cwl-inputs/

export GOOGLE_AUTH_MODE
export PAPI_INPUT_GCS_PREFIX

cromwell::build::run_conformance

cromwell::build::generate_code_coverage
