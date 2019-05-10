#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_docker_environment

# When running locally, an `ADD .` in the Dockerfile might pull in files the local developer does not expect.
# Until there is further review of what goes into the docker image (ex: rendered vault secrets!) do not push it yet.
docker_tag="broadinstitute/cromwell-docker-develop:test-only-do-not-push"

docker build -t "${docker_tag}" scripts/docker-develop

echo "What tests would you like, my dear McMuffins?"

echo "1. Testing for install of sbt"
docker run --rm "${docker_tag}" which sbt

echo "2. Testing sbt assembly"
docker run --rm -v "${PWD}:${PWD}" -w "${PWD}" "${docker_tag}" sbt assembly

echo "3. Testing cloudwell docker compose"

CROMWELL_TAG=develop
export CROMWELL_TAG

docker-compose -f scripts/docker-compose-mysql/docker-compose-cloudwell.yml up --scale cromwell=3 -d

# Give them some time to be ready
sleep 30

# Set the test case
CENTAUR_TEST_FILE=scripts/docker-compose-mysql/test/hello.test
export CENTAUR_TEST_FILE

# Call centaur with our custom test case
sudo sbt "centaur/it:testOnly *ExternalTestCaseSpec"

# Tear everything down
docker-compose -f scripts/docker-compose-mysql/docker-compose-cloudwell.yml down
