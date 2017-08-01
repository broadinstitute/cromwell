#!/usr/bin/env bash

shutdown() {
    cd "${INITIAL_DIR}"
    exit "${EXIT_CODE}"
}

trap "shutdown" EXIT

set -e

EXIT_CODE=1
PROGNAME="$(basename $0)"

usage="
$PROGNAME [-b branch] [-j jar path] [-r rundir] [-c config file] [-t refresh token] [-e excludeTag] [-i testDirPath]

Builds and runs specified branch of Cromwell and runs Centaur against it.

Arguments:
    -b    Branch of Cromwell to test. Mutually exclusive with -j
    -j    Path of a cromwell jar to use. Mutually exclusive with -b
    -r    Directory where script is run (defaults to current directory)
    -c    If supplied, the config file to pass to Cromwell
    -t    Refresh Token that can be passed into the appropriate options file
    -e    If supplied, will exclude tests with this tag
    -p    If supplied, number of tests to be run in parallel. 16 is the default
    -i    If supplied, will run the tests in this directory instead of the standard tests
"

INITIAL_DIR=$(pwd)
RUN_DIR=$(pwd)
TEST_THREAD_COUNT=16

while getopts ":hb:r:c:p:j:t:e:i:" option; do
    case "$option" in
        h) echo "$usage"
            exit
            ;;
        b) CROMWELL_BRANCH="${OPTARG}"
            ;;
        r) RUN_DIR=$OPTARG
            mkdir -p "${RUN_DIR}"
            ;;
        c) CONFIG_STRING="${OPTARG}"
            ;;

        j) CROMWELL_JAR="${OPTARG}"
            ;;
        t) REFRESH_TOKEN=-Dcentaur.optionalToken=$(cat "${OPTARG}")
            ;;
        p) TEST_THREAD_COUNT="${OPTARG}"
            ;;
        i) TEST_CASE_DIR="${OPTARG}"
            ;;
        e) EXCLUDE_TAG+=("${OPTARG}")
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

LOG_DIR="${RUN_DIR}"/logs
ASSEMBLY_LOG=${LOG_DIR}/cromwell_assembly.log
CROMWELL_LOG=${LOG_DIR}/cromwell.log
CENTAUR_LOG=${LOG_DIR}/centaur.log

cd "${RUN_DIR}"
mkdir -p "${LOG_DIR}"

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
    sbt clean assembly >> "${ASSEMBLY_LOG}" 2>&1
    cd ..
    CROMWELL_JAR=$(find "${RUN_DIR}"/cromwell/target/scala-2.* -name "cromwell-*.jar")
fi

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

sbt test:compile
CP=$(sbt "export test:dependency-classpath" --error)

if [ -n "${TEST_CASE_DIR}" ]; then
    RUN_SPECIFIED_TEST_DIR_CMD="-Dcentaur.standardTestCasePath=${TEST_CASE_DIR}"
fi

CENTAUR_CROMWELL_MODE="-Dcentaur.cromwell.mode=jar"
CENTAUR_CROMWELL_JAR="-Dcentaur.cromwell.jar.path=${CROMWELL_JAR}"
CENTAUR_CROMWELL_CONF="-Dcentaur.cromwell.jar.conf=${CONFIG_STRING}"
CENTAUR_CROMWELL_LOG="-Dcentaur.cromwell.jar.log=${CROMWELL_LOG}"
CENTAUR_CROMWELL_RESTART="-Dcentaur.cromwell.jar.withRestart=true"
CENTAUR_CONF="${CENTAUR_CROMWELL_MODE} ${CENTAUR_CROMWELL_JAR} ${CENTAUR_CROMWELL_CONF} ${CENTAUR_CROMWELL_LOG} ${CENTAUR_CROMWELL_RESTART}"

if [[ -n ${EXCLUDE_TAG[*]} ]]; then
    echo "Running Centaur filtering out ${EXCLUDE_TAG[*]} tests"
    EXCLUDE=""
    for val in "${EXCLUDE_TAG[@]}"; do 
        EXCLUDE="-l $val "$EXCLUDE
    done
    TEST_COMMAND="java ${REFRESH_TOKEN} ${RUN_SPECIFIED_TEST_DIR_CMD} ${CENTAUR_CONF} -cp $CP org.scalatest.tools.Runner -R target/scala-2.12/test-classes -oD -PS${TEST_THREAD_COUNT} "$EXCLUDE
else
    echo "Running Centaur with sbt test"
    TEST_COMMAND="java ${REFRESH_TOKEN} ${RUN_SPECIFIED_TEST_DIR_CMD} ${CENTAUR_CONF} -cp $CP org.scalatest.tools.Runner -R target/scala-2.12/test-classes -oD -PS${TEST_THREAD_COUNT}"
fi

eval "${TEST_COMMAND} >> ${CENTAUR_LOG} 2>&1"

if [ $? -eq 0 ]; then
    EXIT_CODE=0
    TEST_STATUS="succeeded"
fi

echo "SBT test $TEST_STATUS, please see logs for more information"
