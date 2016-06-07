#!/bin/bash

shutdown() {
    cd "${INITIAL_DIR}"
    # This will take out the backgrounded Cromwell instance
    pkill -P $$
    exit 0
}

trap "shutdown" EXIT

set -e

PROGNAME="$(basename $0)"

usage="
$PROGNAME [-b branch] [-r rundir] [-c config file] [-p parallelism factor]

Builds and runs specified branch of Cromwell and runs Centaur against it.

Arguments:
    -b    Branch of Cromwell to test
    -r    Directory where script is run (defaults to current directory)
    -c    If supplied, the config file to pass to Cromwell
    -p    Number of simultaneous tests to run. Defaults to 3
"

INITIAL_DIR=$(pwd)

PARALLELISM_FACTOR=3
CROMWELL_BRANCH="develop"
RUN_DIR=$(pwd)
CONFIG_STRING=""

while getopts ":hb:r:c:p:" option; do
    case "$option" in
	h) echo "$usage"
	   exit
	   ;;
	b) CROMWELL_BRANCH=$OPTARG
	   ;;
        r) RUN_DIR=$OPTARG
	   mkdir -p "${RUN_DIR}"
	   ;;
        c) CONFIG_STRING="-Dconfig.file=$OPTARG"
	   ;;
        p) PARALLELISM_FACTOR=$OPTARG
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
shift $((OPTIND - 1))

LOG_DIR=logs
ASSEMBLY_LOG=${LOG_DIR}/cromwell_assembly.log
CROMWELL_LOG=${LOG_DIR}/cromwell.log
CENTAUR_LOG=${LOG_DIR}/centaur.log

cd "${RUN_DIR}"
mkdir -p ${LOG_DIR}

# Build and run cromwell
if [ ! -d "cromwell" ]; then
    echo "Cloning Cromwell"
    git clone https://github.com/broadinstitute/cromwell.git
else
    echo "Using existing Cromwell directory"
fi

cd cromwell
git checkout "${CROMWELL_BRANCH}"
git pull
echo "Building Cromwell"
sbt clean assembly >> ../${ASSEMBLY_LOG} 2>&1

java "${CONFIG_STRING}" -jar target/scala-2.11/cromwell-*.jar server >> ../${CROMWELL_LOG} 2>&1 &

# Build and run centaur
cd ..

if [ ! -d "centaur" ]; then
    echo "Cloning Centaur"
    git clone https://github.com/broadinstitute/centaur.git
else
    echo "Using an existing Centaur directory"
fi

cd centaur

echo "Running Centaur with ${PARALLELISM_FACTOR}-way parallelism"
if ./run_tests_parallel.sh "${PARALLELISM_FACTOR}"  >> ../${CENTAUR_LOG} 2>&1 ; then
    TEST_STATUS="failed"
else
    TEST_STATUS="succeeded"
fi

echo "SBT test $TEST_STATUS, please see logs for more information"
tail -n 5 ../${CENTAUR_LOG}
