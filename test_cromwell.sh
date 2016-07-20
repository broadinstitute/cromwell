#!/bin/bash

shutdown() {
    cd "${INITIAL_DIR}"
    # This will take out the backgrounded Cromwell instance
    pkill -P $$
    exit "${EXIT_CODE}"
}

trap "shutdown" EXIT

set -e

EXIT_CODE=1
PROGNAME="$(basename $0)"

usage="
$PROGNAME [-b branch] [-j jar path] [-r rundir] [-c config file] [-p parallelism factor] [-t refresh token]

Builds and runs specified branch of Cromwell and runs Centaur against it.

Arguments:
    -b    Branch of Cromwell to test. Mutually exclusive with -j
    -j    Path of a cromwell jar to use. Mutually exclusive with -b
    -r    Directory where script is run (defaults to current directory)
    -c    If supplied, the config file to pass to Cromwell
    -p    Number of simultaneous tests to run. Defaults to 3
    -t    Refresh Token that can be passed into the appropriate options file
"

INITIAL_DIR=$(pwd)
RUN_DIR=$(pwd)

while getopts ":hb:r:c:p:j:t:" option; do
    case "$option" in
        h) echo "$usage"
            exit
            ;;
        b) CROMWELL_BRANCH="${OPTARG}"
            ;;
        r) RUN_DIR=$OPTARG
            mkdir -p "${RUN_DIR}"
            ;;
        c) CONFIG_STRING=-Dconfig.file="${OPTARG}"
            ;;
        p) PARALLELISM_FACTOR="${OPTARG}"
            ;;
        j) CROMWELL_JAR="${OPTARG}"
            ;;
        t) REFRESH_TOKEN=-Dcentaur.optionalToken="${OPTARG}"
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
shift "$((OPTIND - 1))"

LOG_DIR=logs
ASSEMBLY_LOG=${LOG_DIR}/cromwell_assembly.log
CROMWELL_LOG=${LOG_DIR}/cromwell.log
CENTAUR_LOG=${LOG_DIR}/centaur.log

cd "${RUN_DIR}"
mkdir -p ${LOG_DIR}

if [[ -n ${CROMWELL_BRANCH} ]]; then
    if [[ -n ${CROMWELL_JAR} ]]; then
        echo "Do not specify both a branch and a jar" >&2
        exit 1
    fi

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
    cd ..
    CROMWELL_JAR="cromwell/target/scala-2.11/cromwell-*.jar"
fi

echo "Starting Cromwell, jar is ${CROMWELL_JAR}"
java "${CONFIG_STRING}" -jar "${CROMWELL_JAR}" server >> "${CROMWELL_LOG}" 2>&1 &

# Build and run centaur
cd "${RUN_DIR}"

if [[ ! -d "centaur" ]]; then
    echo "Cloning Centaur"
    git clone https://github.com/broadinstitute/centaur.git
else
    echo "Using an existing Centaur directory"
fi

cd centaur
TEST_STATUS="failed"

if [[ -n ${PARALLELISM_FACTOR} ]]; then
    echo "Running Centaur with ${PARALLELISM_FACTOR}-way parallelism"
    TEST_COMMAND="./run_tests_parallel.sh ${PARALLELISM_FACTOR}"
else
    echo "Running Centaur with sbt test"
    echo "About to run centaur as sbt ${REFRESH_TOKEN} test"
    TEST_COMMAND="sbt ${REFRESH_TOKEN} test"
fi

${TEST_COMMAND}  >> ../${CENTAUR_LOG} 2>&1

if [ $? -eq 0 ]; then
    EXIT_CODE=0
    TEST_STATUS="succeeded"
fi

echo "SBT test $TEST_STATUS, please see logs for more information"
tail -n 5 ../${CENTAUR_LOG}
