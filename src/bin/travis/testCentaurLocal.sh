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

sbt assembly
CROMWELL_JAR=$(find "$(pwd)/target/scala-2.11" -name "cromwell-*.jar")
LOCAL_CONF="$(pwd)/src/bin/travis/resources/local_centaur.conf"
git clone https://github.com/broadinstitute/centaur.git
cd centaur
git checkout ${CENTAUR_BRANCH}
cd ..
# All tests use ubuntu:latest - make sure it's there before starting the tests 
# because pulling the image during some of the tests would cause them to fail 
# (specifically output_redirection which expects a specific value in stderr)
docker pull ubuntu:latest
centaur/test_cromwell.sh -j"${CROMWELL_JAR}" -c${LOCAL_CONF}
