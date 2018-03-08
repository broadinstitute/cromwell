#!/usr/bin/env bash

if [ "$TRAVIS_SECURE_ENV_VARS" = "false" ]; then
    echo "************************************************************************************************"
    echo "************************************************************************************************"
    echo "**                                                                                            **"
    echo "**  WARNING: Encrypted keys are unavailable to automatically test JES with centaur. Exiting.  **"
    echo "**                                                                                            **"
    echo "************************************************************************************************"
    echo "************************************************************************************************"
    exit 0
fi

printTravisHeartbeat() {
    # Sleep one minute between printouts, but don't zombie for more than two hours
    for ((i=0; i < 180; i++)); do
        sleep 60
        printf "…"
    done &
    TRAVIS_HEARTBEAT_PID=$!
}

cromwellLogTail() {
 (
   while [ ! -f logs/cromwell.log ];
   do
     sleep 2
     printf "(Cr)"
   done
   tail -f logs/cromwell.log &
   CROMWELL_LOG_TAIL_PID=$!
 ) &
 CROMWELL_LOG_WAIT_PID=$!
}

centaurLogTail() {
 (
   while [ ! -f logs/centaur.log ];
   do
     sleep 2
     printf "(Ce)"
   done
   tail -f logs/centaur.log &
   CENTAUR_LOG_TAIL_PID=$!
 ) &
 CENTAUR_LOG_WAIT_PID=$!
}

killTravisHeartbeat() {
    if [ -n "${TRAVIS_HEARTBEAT_PID+set}" ]; then
        kill ${TRAVIS_HEARTBEAT_PID} || true
    fi
}

killCromwellLogTail() {
    if [ -n "${CROMWELL_LOG_TAIL_PID+set}" ]; then
        kill ${CROMWELL_LOG_TAIL_PID} || true
    else
        if [ -n "${CROMWELL_LOG_WAIT_PID+set}" ]; then
            kill ${CROMWELL_LOG_WAIT_PID} || true
        fi
    fi
}

killCentaurLogTail() {
    if [ -n "${CENTAUR_LOG_TAIL_PID+set}" ]; then
        kill ${CENTAUR_LOG_TAIL_PID} || true
    else
        if [ -n "${CENTAUR_LOG_WAIT_PID+set}" ]; then
            kill ${CENTAUR_LOG_WAIT_PID} || true
        fi
    fi
}

exitScript() {
    echo "CENTAUR LOG"
    cat logs/centaur.log
    killTravisHeartbeat
    killCromwellLogTail
    killCentaurLogTail
}

trap exitScript EXIT
trap exitScript TERM
cromwellLogTail
centaurLogTail
printTravisHeartbeat

set -x
set -e

PROGNAME="$(basename "$0")"
RUN_INTEGRATION_TESTS=0

usage="
$PROGNAME [-i ]

Builds and runs Cromwell and runs Centaur against it.

Arguments:
    -i    Flag that if supplied, will run centaur integration tests instead of standardtests
"

while getopts ":hi" option; do
    case "$option" in
        h) echo "$usage"
            exit
            ;;
        i) RUN_INTEGRATION_TESTS=1
            ;;
        :) printf "Missing argument for -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1
            ;;
        \?) printf "Illegal option: -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1
            ;;
        esac
done

# TURN OFF LOGGING WHILE WE TALK TO DOCKER/VAULT
set +x

# Login to docker to access the dsde-toolbox
docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"

# Login to vault to access secrets
docker run --rm \
    -v $HOME:/root:rw \
    broadinstitute/dsde-toolbox \
    vault auth "$JES_TOKEN" < /dev/null > /dev/null && echo vault auth success

set -x

# Render secrets
docker run --rm \
    -v $HOME:/root:rw \
    -v $PWD/src/bin/travis/resources:/working \
    -v $PWD:/output \
    -e ENVIRONMENT=not_used \
    -e INPUT_PATH=/working \
    -e OUT_PATH=/output \
    broadinstitute/dsde-toolbox render-templates.sh

ASSEMBLY_LOG_LEVEL=error ENABLE_COVERAGE=true sbt assembly --error
CROMWELL_JAR=$(find "$(pwd)/server/target/scala-2.12" -name "cromwell-*.jar")
JES_CONF="$(pwd)/jes_centaur.conf"
GOOGLE_AUTH_MODE="service-account"
GOOGLE_REFRESH_TOKEN_PATH="$(pwd)/jes_refresh_token.txt"
GOOGLE_SERVICE_ACCOUNT_JSON="$(pwd)/cromwell-service-account.json"

# pass integration directory to the inputs json otherwise remove it from the inputs file
if [ $RUN_INTEGRATION_TESTS -ne 1 ]; then
    INTEGRATION_TESTS=""
else
    INTEGRATION_TESTS="-i$INTEGRATION_TESTS_DIR"
fi

# All tests use ubuntu:latest - make sure it's there before starting the tests
# because pulling the image during some of the tests would cause them to fail
# (specifically output_redirection which expects a specific value in stderr)
docker pull ubuntu:latest

# Export variables used in conf files
export GOOGLE_AUTH_MODE
export GOOGLE_REFRESH_TOKEN_PATH
export GOOGLE_SERVICE_ACCOUNT_JSON

centaur/test_cromwell.sh \
  -j${CROMWELL_JAR} \
  -g \
  -c${JES_CONF} \
  -elocaldockertest \
  -p100 \
  $INTEGRATION_TESTS

if [ "$TRAVIS_EVENT_TYPE" != "cron" ]; then
    sbt coverageReport --warn
    sbt coverageAggregate --warn
    bash <(curl -s https://codecov.io/bash) >/dev/null
fi
