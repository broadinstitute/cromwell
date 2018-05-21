#!/usr/bin/env bash

set -e
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_everyone_environment

ENABLE_COVERAGE=true sbt \
    -Dbackend.providers.Local.config.filesystems.local.localization.0=copy \
    +clean +nointegration:test

cromwell::build::generate_code_coverage

cromwell::build::publish_artifacts
