#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

export CROMWELL_BUILD_REQUIRES_SECURE=true

# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

docker run --rm \
           -v "${CROMWELL_BUILD_ROOT_DIRECTORY}/CromwellRefdiskManifestCreator:/CromwellRefdiskManifestCreator" \
           maven:3.6.3-openjdk-11 /bin/bash -c "
cd /CromwellRefdiskManifestCreator
mvn test
"
