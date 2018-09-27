#!/usr/bin/env bash

set -e
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

case "${CROMWELL_BUILD_PROVIDER}" in
    "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
        CROMWELL_SBT_TEST_EXCLUDE_TAGS="AwsTest,CromwellIntegrationTest,GcsIntegrationTest"
        ;;
    "${CROMWELL_BUILD_PROVIDER_JENKINS}")
        CROMWELL_SBT_TEST_EXCLUDE_TAGS="AwsTest,CromwellIntegrationTest,DockerTest,GcsIntegrationTest"
        ;;
    *)
        # Use the full list of excludes listed in Testing.scala
        CROMWELL_SBT_TEST_EXCLUDE_TAGS=""
        ;;
esac
export CROMWELL_SBT_TEST_EXCLUDE_TAGS

sbt -Dbackend.providers.Local.config.filesystems.local.localization.0=copy coverage test

cromwell::build::generate_code_coverage

cromwell::build::publish_artifacts
