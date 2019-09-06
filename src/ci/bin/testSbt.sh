#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR=1

case "${CROMWELL_BUILD_PROVIDER}" in
    "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
        CROMWELL_SBT_TEST_EXCLUDE_TAGS="AwsTest,CromwellIntegrationTest,GcsIntegrationTest"
        CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR=2
        ;;
    "${CROMWELL_BUILD_PROVIDER_JENKINS}")
        CROMWELL_SBT_TEST_EXCLUDE_TAGS="AwsTest,CromwellIntegrationTest,DockerTest,GcsIntegrationTest"
        CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR=10
        ;;
    *)
        # Use the full list of excludes listed in Testing.scala
        CROMWELL_SBT_TEST_EXCLUDE_TAGS=""
        ;;
esac
export CROMWELL_SBT_TEST_EXCLUDE_TAGS
export CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR

echo "Installed Java version:"
java -version

echo "Contents of /etc/sbt/jvmopts:"
cat /etc/sbt/jvmopts
echo "Contents of /etc/sbt/sbtopts"
cat /etc/sbt/sbtopts

echo "Running Java memory experiment..."
javac TotalMemory.java
java -Xms2048M -Xmx2048M -Xss6M -XX:MaxPermSize=512M -Xmx4g -XX:MaxMetaspaceSize=2g -Xss8m -cp . TotalMemory

echo "Starting SBT..."

# `-v` prints SBT's environment variables and then carries on executing the requested command as usual
sbt -v -Dsbt.classloader.close=false -Dakka.test.timefactor=${CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR} -Dbackend.providers.Local.config.filesystems.local.localization.0=copy coverage test

cromwell::build::generate_code_coverage

cromwell::build::publish_artifacts
