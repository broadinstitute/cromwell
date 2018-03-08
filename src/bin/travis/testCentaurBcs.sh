#!/usr/bin/env bash

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

ASSEMBLY_LOG_LEVEL=error ENABLE_COVERAGE=true sbt assembly --error
CROMWELL_JAR=$(find "$(pwd)/server/target/scala-2.12" -name "cromwell-*.jar")
BCS_CONF="$(pwd)/src/bin/travis/resources//bcs_centaur.conf"

# All tests use ubuntu:latest - make sure it's there before starting the tests
# because pulling the image during some of the tests would cause them to fail 
# (specifically output_redirection which expects a specific value in stderr)
docker pull ubuntu:latest

# We exclude these cases because that BCS backend can not support call caching, 
# docker image in docker hub, glob files and cwl currently
centaur/test_cromwell.sh -j"${CROMWELL_JAR}" -g -c${BCS_CONF} \
-e call_cache_capoeira_local \
-e non_root_default_user \
-e non_root_specified_user \
-e docker_hash_dockerhub \
-e valid_labels \
-e scattergather \
-e globbingscatter \
-e hello \
-e three_step__subwf_cwl \
-e invalidate_bad_caches_local \
-e three_step_cwl \
-e hello_yaml \
-e inline_file \
-e docker_hash_quay \
-e globbingindex \
-e dontglobinputs \
-e space \
-e docker_hash_gcr \
-e curl \
-e globbingbehavior \
-e multiline_command_line \
-e dont_cache_to_failed_jobs \
-e readfromcachefalse \
-e floating_tags \
-e cachewithinwf \
-e cachebetweenwf \
-e writetocache \
-e test_file_outputs_from_input \
-e array_io \
-e bad_file_string \
-e sub_workflow_interactions \
-e lots_of_inputs \
-e no_new_calls \
-e bad_output_task 

if [ "$TRAVIS_EVENT_TYPE" != "cron" ]; then
    sbt coverageReport --warn
    sbt coverageAggregate --warn
    bash <(curl -s https://codecov.io/bash) >/dev/null
fi
