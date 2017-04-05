#!/usr/bin/env bash

printTravisHeartbeat() {
    # Sleep one minute between printouts, but don't zombie for more than two hours
    for ((i=0; i < 120; i++)); do
        sleep 60
        printf "…"
    done &
    TRAVIS_HEARTBEAT_PID=$!
}

killTravisHeartbeat() {
    if [ -n "${TRAVIS_HEARTBEAT_PID+set}" ]; then
        kill ${TRAVIS_HEARTBEAT_PID} || true
    fi
}

exitScript() {
    echo "TES LOG"
    cat logs/tes.log
    echo "CROMWELL LOG"
    cat logs/cromwell.log
    echo "CENTAUR LOG"
    cat logs/centaur.log
    killTravisHeartbeat
}

trap exitScript EXIT
printTravisHeartbeat

set -x
set -e

sbt assembly
CROMWELL_JAR=$(find "$(pwd)/target/scala-2.11" -name "cromwell-*.jar")
TES_CENTAUR_CONF="$(pwd)/src/bin/travis/resources/tes_centaur.conf"
git clone https://github.com/broadinstitute/centaur.git
cd centaur
git checkout ${CENTAUR_BRANCH}
cd ..

TES_CONF="$(pwd)/src/bin/travis/resources/tes.conf"
git clone https://github.com/ohsu-comp-bio/funnel.git
cd funnel
git checkout 7adbf0b
make
cd ..
mkdir logs
nohup funnel/bin/tes-server -config ${TES_CONF} > logs/tes.log 2>&1 &


# All tests use ubuntu:latest - make sure it's there before starting the tests 
# because pulling the image during some of the tests would cause them to fail 
# (specifically output_redirection which expects a specific value in stderr)
docker pull ubuntu:latest
centaur/test_cromwell.sh -j"${CROMWELL_JAR}" -c${TES_CENTAUR_CONF} -elocaldockertest
