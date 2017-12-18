#!/usr/bin/env bash

set -e

sudo -H pip install --upgrade pip
sudo -H pip install cwltest

# -- BEGIN PAPI related conf generation

# TURN OFF LOGGING WHILE WE TALK TO DOCKER/VAULT
set +x

# Login to docker to access the dsde-toolbox
docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"

# Login to vault to access secrets
docker run --rm \
    -v "$HOME:/root:rw" \
    broadinstitute/dsde-toolbox \
    vault auth "$JES_TOKEN" < /dev/null > /dev/null && echo vault auth success

set -x

# Render secrets
docker run --rm \
    -v "$HOME:/root:rw" \
    -v "$PWD/src/bin/travis/resources:/working" \
    -v "$PWD:/output" \
    -e ENVIRONMENT=not_used \
    -e INPUT_PATH=/working \
    -e OUT_PATH=/output \
    broadinstitute/dsde-toolbox render-templates.sh

JES_CONF="$(pwd)/jes_centaur.conf"

# -- END PAPI related conf generation

# Turn off call caching as hashing doesn't work since it sees local and not GCS paths.
sed -i '/^call-caching\s*/{N;s/enabled.*/  enabled: false/;}' ${JES_CONF}

ENABLE_COVERAGE=true sbt assembly
CROMWELL_JAR=$(find "$(pwd)/target/scala-2.12" -name "cromwell-*.jar")
CENTAUR_CWL_RUNNER="$(pwd)/centaurCwlRunner/src/bin/centaur-cwl-runner.bash"

git clone --depth 1 https://github.com/common-workflow-language/common-workflow-language.git

CONFORMANCE_EXPECTED_FAILURES=$(pwd)/src/bin/travis/resources/conformance_expected_failures.txt

shutdownCromwell() {
    if [ -z "${CROMWELL_PID}" ]; then
        kill "${CROMWELL_PID}"
    fi
}

trap "shutdownCromwell" EXIT
CURRENT_DIR=$(pwd)

java \
  -Dconfig.file="$JES_CONF" \
  -Dpapi.default-input-gcs-prefix=gs://centaur-cwl-conformance/cwl-inputs/ \
  -Dcwl.default-docker-image="ubuntu:latest" \
  -Dsystem.new-workflow-poll-rate=1 \
  -jar "${CROMWELL_JAR}" server &

CROMWELL_PID=$$

sleep 30

UNEXPECTED_PASS=()
UNEXPECTED_FAIL=()

cd common-workflow-language
TEST_COUNT=$(./run_test.sh RUNNER="${CENTAUR_CWL_RUNNER}" -l | grep -c '^\[')

echo Test count is "${TEST_COUNT}"

TEST_PASSING=0

set +e
# Only run tests that are expected to succeed as PAPI CWL conformance is not quick.
for TEST_NUMBER in $(perl -ne '/#(\d+)/ && print "$1 "' "${CONFORMANCE_EXPECTED_FAILURES}"); do
    # Check if test is supposed to fail
    grep -q '^'"${TEST_NUMBER}"'$' "${CONFORMANCE_EXPECTED_FAILURES}"
    TEST_IN_EXPECTED_FAILED=$?

    # Run the test
    ./run_test.sh RUNNER="${CENTAUR_CWL_RUNNER}" -n"${TEST_NUMBER}"
    TEST_RESULT_CODE=$?

    if [ ${TEST_RESULT_CODE} -eq 0 ]; then
        TEST_PASSING=$((TEST_PASSING + 1))
    fi

    # Check for unexpected results
    if [ ${TEST_IN_EXPECTED_FAILED} -eq 0 ] && [ ${TEST_RESULT_CODE} -eq 0 ]; then
        UNEXPECTED_PASS+=("${TEST_NUMBER}")
    elif [ ! ${TEST_IN_EXPECTED_FAILED} -eq 0 ] && [ ! ${TEST_RESULT_CODE} -eq 0 ]; then
        UNEXPECTED_FAIL+=("${TEST_NUMBER}")
    fi
done
set -e

cd "$CURRENT_DIR"

echo Conformance percentage at "$(echo 100 '*' ${TEST_PASSING} '/' "${TEST_COUNT}" | bc)%"

if [ ! $(( ${#UNEXPECTED_PASS[@]} + ${#UNEXPECTED_FAIL[@]} )) -eq 0 ]; then
    printf 'Unexpected passing tests: (%s)\n' "${UNEXPECTED_PASS[*]}"
    printf 'Unexpected failing tests: (%s)\n' "${UNEXPECTED_FAIL[*]}"
    printf 'Does %s need to be updated?\n' "${CONFORMANCE_EXPECTED_FAILURES}"
    exit 1
fi

if [ "$TRAVIS_EVENT_TYPE" != "cron" ]; then
    sbt coverageReport --warn
    sbt coverageAggregate --warn
    bash <(curl -s https://codecov.io/bash) >/dev/null
fi
