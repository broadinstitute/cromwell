#!/usr/bin/env bash

if [ "$TRAVIS_SECURE_ENV_VARS" = "false" ]; then
    echo "************************************************************************************************"
    echo "************************************************************************************************"
    echo "**                                                                                            **"
    echo "**  WARNING: Encrypted keys are unavailable to automatically test BCS with centaur. Exiting.  **"
    echo "**                                                                                            **"
    echo "************************************************************************************************"
    echo "************************************************************************************************"
    exit 0
fi

printTravisHeartbeat() {
    # Sleep one minute between printouts, but don't zombie for more than two hours
    for ((i=0; i < 120; i++)); do
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
BCS_CONF="$(pwd)/bcs_centaur.conf"

# All tests that run on the hardwired Local backend use ubuntu:latest - make sure it's there before starting the tests
# because pulling the image during some of the tests would cause them to fail 
# (specifically output_redirection which expects a specific value in stderr)
docker pull ubuntu:latest

# https://github.com/broadinstitute/cromwell/issues/3522
# https://github.com/broadinstitute/cromwell/issues/3523
# https://github.com/broadinstitute/cromwell/issues/3524
exclude_known_bugs=" \
  -e bad_file_string \
  -e bad_output_task \
  -e tmp_dir \
"

# https://github.com/broadinstitute/cromwell/issues/3518
exclude_docker_tests=" \
  -e curl \
  -e docker_hash_dockerhub \
  -e docker_hash_gcr \
  -e docker_hash_quay \
  -e dont_cache_to_failed_jobs \
  -e hello \
  -e hello_yaml \
  -e inline_file \
  -e inline_file_custom_entryname \
  -e iwdr_input_string \
  -e iwdr_input_string_function \
  -e non_root_default_user \
  -e three_step_cwl \
"

# https://github.com/broadinstitute/cromwell/issues/3519
exclude_glob_tests=" \
  -e cwl_glob_sort \
  -e cwl_interpolated_strings \
  -e dontglobinputs \
  -e globbingbehavior \
  -e globbingindex \
  -e globbingscatter \
  -e lots_of_inputs \
  -e space \
  -e wdl_empty_glob \
"

centaur/test_cromwell.sh \
  -j"${CROMWELL_JAR}" \
  -g \
  -c${BCS_CONF} \
  -elocaldockertest \
  -p100 \
  -t1m \
  $exclude_known_bugs \
  $exclude_docker_tests \
  $exclude_glob_tests \

if [ "$TRAVIS_EVENT_TYPE" != "cron" ]; then
    sbt coverageReport --warn
    sbt coverageAggregate --warn
    bash <(curl -s https://codecov.io/bash) >/dev/null
fi
