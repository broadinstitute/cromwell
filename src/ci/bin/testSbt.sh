#!/usr/bin/env bash

set -e
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

sbt \
    -Dbackend.providers.Local.config.filesystems.local.localization.0=copy \
    coverage nointegration:test

cromwell::build::generate_code_coverage

cromwell::build::publish_artifacts
