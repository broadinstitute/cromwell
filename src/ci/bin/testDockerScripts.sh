#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

# When running locally, an `ADD .` in the Dockerfile might pull in files the local developer does not expect.
# Until there is further review of what goes into the docker image (ex: rendered vault secrets!) do not push it yet.
docker_tag="broadinstitute/cromwell-docker-develop:test-only-do-not-push"

docker build -t "${docker_tag}" scripts/docker-develop

echo "What tests would you like, my dear McMuffins?"

echo "1. Testing for install of sbt"
docker run --rm "${docker_tag}" which sbt

echo "2. Testing sbt assembly"
docker run \
    --rm \
    --volume "${CROMWELL_BUILD_ROOT_DIRECTORY}:${CROMWELL_BUILD_ROOT_DIRECTORY}" \
    --workdir "${CROMWELL_BUILD_ROOT_DIRECTORY}" \
    "${docker_tag}" \
    sbt --warn assembly
