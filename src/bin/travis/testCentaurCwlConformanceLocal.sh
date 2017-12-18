#!/usr/bin/env bash

set -e

sudo -H pip install --upgrade pip
sudo -H pip install cwltest

ENABLE_COVERAGE=true sbt assembly
CROMWELL_JAR=$(find "$(pwd)/target/scala-2.12" -name "cromwell-*.jar")
CENTAUR_CWL_RUNNER="$(pwd)/centaurCwlRunner/src/bin/centaur-cwl-runner.bash"

git clone --depth 1 https://github.com/common-workflow-language/common-workflow-language.git
CWL_TEST_DIR=$(pwd)/common-workflow-language/v1.0/v1.0
CONFORMANCE_EXPECTED_FAILURES=$(pwd)/src/bin/travis/resources/conformance_expected_failures.txt

shutdownCromwell() {
    if [ -z "${CROMWELL_PID}" ]; then
        kill ${CROMWELL_PID}
    fi
}

trap "shutdownCromwell" EXIT

# Start the Cromwell server in the directory containing input files so it can access them via their relative path
CURRENT_DIR=$(pwd)
cd $CWL_TEST_DIR
java -Dsystem.new-workflow-poll-rate=1 -Dbackend.providers.Local.config.script-epilogue="sleep 5 && sync" -jar ${CROMWELL_JAR} server &
CROMWELL_PID=$$
cd $CURRENT_DIR

sleep 5

UNEXPECTED_PASS=()
UNEXPECTED_FAIL=()

cd common-workflow-language
TEST_COUNT=$(./run_test.sh RUNNER="${CENTAUR_CWL_RUNNER}" -l | grep '^\[' | wc -l)

echo Test count is ${TEST_COUNT}

TEST_PASSING=0

set +e
for TEST_NUMBER in $(seq ${TEST_COUNT}); do
    # Check if test is supposed to fail
    grep -q '^'${TEST_NUMBER}'$' ${CONFORMANCE_EXPECTED_FAILURES}
    TEST_IN_EXPECTED_FAILED=$?

    # Run the test
    ./run_test.sh RUNNER="${CENTAUR_CWL_RUNNER}" -n${TEST_NUMBER}
    TEST_RESULT_CODE=$?

    if [ ${TEST_RESULT_CODE} -eq 0 ]; then
        TEST_PASSING=$(($TEST_PASSING + 1))
    fi

    # Check for unexpected results
    if [ ${TEST_IN_EXPECTED_FAILED} -eq 0 ] && [ ${TEST_RESULT_CODE} -eq 0 ]; then
        UNEXPECTED_PASS+=(${TEST_NUMBER})
    elif [ ! ${TEST_IN_EXPECTED_FAILED} -eq 0 ] && [ ! ${TEST_RESULT_CODE} -eq 0 ]; then
        UNEXPECTED_FAIL+=(${TEST_NUMBER})
    fi
done
set -e

cd $CURRENT_DIR

echo Conformance percentage at "$(echo 100 '*' ${TEST_PASSING} '/' ${TEST_COUNT} | bc)%"

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
