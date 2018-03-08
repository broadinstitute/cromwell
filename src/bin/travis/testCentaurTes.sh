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

TES_CENTAUR_CONF="$(pwd)/src/bin/travis/resources/tes_centaur.conf"
FUNNEL_CONF="$(pwd)/src/bin/travis/resources/funnel.conf"

wget https://github.com/ohsu-comp-bio/funnel/releases/download/0.5.0/funnel-linux-amd64-0.5.0.tar.gz
tar xzf funnel-linux-amd64-0.5.0.tar.gz
FUNNEL_PATH="$(pwd)/funnel"

mkdir logs
nohup "${FUNNEL_PATH}" server run --config "${FUNNEL_CONF}" > logs/funnel.log 2>&1 &

# All tests use ubuntu:latest - make sure it's there before starting the tests
# because pulling the image during some of the tests would cause them to fail 
# (specifically output_redirection which expects a specific value in stderr)
docker pull ubuntu:latest

# The following tests are skipped:
#
# call_cache_capoeira_local: fails on task 'read_files_without_docker' since the 'docker' runtime key is required for this backend
# lots_of_inputs:            Funnel mounts in each input separately, this task surpasses the docker limit for volumes
# no_new_calls:              TES does not support checking job status after restart, and cannot tell if shouldSucceed is done or failed
# non_root_specified_user:   TES doesn't support switching users in the image
# write_lines_files:         all inputs are read-only in TES

centaur/test_cromwell.sh \
-j ${CROMWELL_JAR} \
-g \
-c ${TES_CENTAUR_CONF} \
-e call_cache_capoeira_local \
-e lots_of_inputs \
-e no_new_calls \
-e non_root_default_user \
-e non_root_specified_user \
-e write_lines_files \

if [ "$TRAVIS_EVENT_TYPE" != "cron" ]; then
    sbt coverageReport --warn
    sbt coverageAggregate --warn
    bash <(curl -s https://codecov.io/bash) >/dev/null
fi
