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
    echo "FUNNEL LOG"
    cat logs/funnel.log
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

sudo apt-get update -qq
sudo apt-get install -qq mysql-server-5.6 mysql-client-5.6 mysql-client-core-5.6
docker pull ubuntu:latest
mysql -u root -e "SET GLOBAL sql_mode = 'STRICT_ALL_TABLES';"
mysql -u root -e "CREATE DATABASE IF NOT EXISTS cromwell_test;"
mysql -u root -e "CREATE USER 'travis'@'localhost' IDENTIFIED BY '';"
mysql -u root -e "GRANT ALL PRIVILEGES ON cromwell_test . * TO 'travis'@'localhost';"

WORKDIR=$(pwd)

sbt assembly
CROMWELL_JAR=$(find "$(pwd)/target/scala-2.11" -name "cromwell-*.jar")
TES_CENTAUR_CONF="$(pwd)/src/bin/travis/resources/tes_centaur.conf"
git clone https://github.com/broadinstitute/centaur.git
cd centaur
git checkout ${CENTAUR_BRANCH}
cd $WORKDIR


FUNNEL_CONF="$(pwd)/src/bin/travis/resources/funnel.conf"
wget https://storage.googleapis.com/golang/go1.8.1.linux-amd64.tar.gz
tar xfz go1.8.1.linux-amd64.tar.gz
export GOROOT=$WORKDIR/go
mkdir go-lib
export GOPATH=$WORKDIR/go-lib
go get github.com/ohsu-comp-bio/funnel
cd $GOPATH/src/github.com/ohsu-comp-bio/funnel
git checkout 04c5e03
make
cd $WORKDIR
mkdir logs
nohup $GOPATH/bin/funnel server --config ${FUNNEL_CONF} > logs/funnel.log 2>&1 &


# All tests use ubuntu:latest - make sure it's there before starting the tests 
# because pulling the image during some of the tests would cause them to fail 
# (specifically output_redirection which expects a specific value in stderr)
docker pull ubuntu:latest

# The following tests are skipped:
#
# non_root_specified_user:   TES doesn't support switching users in the image
# write_lines_files:         all inputs are read-only in TES
# lots_of_inputs:            Funnel mounts in each input separately, this task surpasses the docker limit for volumes
# call_cache_capoeira_local: fails on task 'read_files_without_docker' since the 'docker' runtime key is required for this backend
#
centaur/test_cromwell.sh \
-j ${CROMWELL_JAR} \
-c ${TES_CENTAUR_CONF} \
-e non_root_specified_user \
-e write_lines_files \
-e lots_of_inputs \
-e call_cache_capoeira_local
