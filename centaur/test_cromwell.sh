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
$PROGNAME [-h] [-b branch] [-j jar path] [-g] [-r rundir] [-l logdir] [-c cromwell config file] [-n centaur config file] [-t refresh token file] [-s service account json] [-i includeTag] [-e excludeTag] [-d testDirPath]

Builds and runs specified branch of Cromwell and runs Centaur against it.

Arguments:
    -h    Displays this help message and exits
    -b    Branch of Cromwell to test. Mutually exclusive with -j
    -j    Path of a Cromwell jar to use. Mutually exclusive with -b
    -r    Directory where script is run (defaults to current directory)
    -l    Directory where logs are written (defaults to logs under the current directory)
    -g    Generate code coverage output for the Centaur main classes
    -c    If supplied, the config file to pass to Cromwell
    -n    If supplied, the config file to pass to Centaur
    -t    If supplied, the timeout for request-plus-response from Centaur to Cromwell
    -i    If supplied, will include tests with this tag
    -e    If supplied, will exclude tests with this tag
    -s    If supplied, will run only the specified suite
    -p    If supplied, number of tests to be run in parallel. 16 is the default
    -d    If supplied, will run the tests in this directory instead of the standard tests
"

INITIAL_DIR=$(pwd)
RUN_DIR=$(pwd)
LOG_DIR="${RUN_DIR}"/logs
TEST_THREAD_COUNT=16
CENTAUR_SBT_COVERAGE=false
CROMWELL_TIMEOUT=10s
SUITE=""

while getopts ":hb:r:l:c:n:p:j:gt:i:e:s:d:" option; do
    case "$option" in
        h) echo "$usage"
            exit
            ;;
        r) RUN_DIR=$OPTARG
            ;;
        l) LOG_DIR=$OPTARG
            ;;
        n) CENTAUR_CONFIG_STRING="${OPTARG}"
            ;;
        g) CENTAUR_SBT_COVERAGE=true
            ;;
        b) CROMWELL_BRANCH="${OPTARG}"
            ;;
        j) CROMWELL_JAR="${OPTARG}"
            ;;
        c) CROMWELL_CONFIG_STRING="${OPTARG}"
            ;;
        t) CROMWELL_TIMEOUT="${OPTARG}"
            ;;
        p) TEST_THREAD_COUNT="${OPTARG}"
            ;;
        d) TEST_CASE_DIR="${OPTARG}"
            ;;
        i) INCLUDE_TAG+=("${OPTARG}")
            ;;
        e) EXCLUDE_TAG+=("${OPTARG}")
            ;;
        s) SUITE="${OPTARG}"
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

ASSEMBLY_LOG=${LOG_DIR}/cromwell_assembly.log
CROMWELL_LOG=${LOG_DIR}/cromwell.log
CENTAUR_LOG=${LOG_DIR}/centaur.log

mkdir -p "${RUN_DIR}"
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
    echo "Checking out CROMWELL_BRANCH $CROMWELL_BRANCH"
    git checkout "${CROMWELL_BRANCH}"
    git pull
    echo "Building Cromwell"
    sbt assembly >> "${ASSEMBLY_LOG}" 2>&1
    cd ..
    # This is the "branch" logic but sets the CROMWELL_JAR to be used in either the "branch" or "jar" use cases.
    # Note that this may not be necessary in the docker-compose use case.
    CROMWELL_JAR=$(find "${RUN_DIR}"/cromwell/target/scala-2.* -name "cromwell-*.jar")
fi

# Build and run centaur
cd "${RUN_DIR}"

TEST_STATUS="failed"

if [[ "${CENTAUR_SBT_COVERAGE}" == "true" ]]; then
    sbt coverage centaur/it:compile
    CP=$(sbt coverage "export centaur/it:dependencyClasspath" -error)
else
    sbt centaur/it:compile
    CP=$(sbt "export centaur/it:dependencyClasspath" -error)
fi

# Add the it-classes folder to the classpath to ensure logback configuration files are picked up.
CP="${CP}:${RUN_DIR}/centaur/target/scala-2.12/it-classes"

# This is set in cromwell::private::create_centaur_variables
if [ -n "${CENTAUR_CONFIG_STRING}" ]; then
    CENTAUR_CONF="-Dconfig.file=${CENTAUR_CONFIG_STRING}"
else
    if [ -n "${TEST_CASE_DIR}" ]; then
        RUN_SPECIFIED_TEST_DIR_CMD="-Dcentaur.standardTestCasePath=${TEST_CASE_DIR}"
    fi

    # Regardless of whether this script was invoked with the "branch" or "jar" option, Centaur is always run in "jar" mode.
    # "branch" or "jar" only controls whether this script was handed a jar or had to build one above.
    CENTAUR_CROMWELL_MODE="-Dcentaur.cromwell.mode=jar"
    CENTAUR_CROMWELL_JAR="-Dcentaur.cromwell.jar.path=${CROMWELL_JAR}"
    CENTAUR_CROMWELL_CONF="-Dcentaur.cromwell.jar.conf=${CROMWELL_CONFIG_STRING}"
    CENTAUR_CROMWELL_LOG="-Dcentaur.cromwell.jar.log=${CROMWELL_LOG}"
    CENTAUR_CROMWELL_RESTART="-Dcentaur.cromwell.jar.withRestart=true"
    CENTAUR_SEND_RECEIVE_TIMEOUT="-Dcentaur.sendReceiveTimeout='${CROMWELL_TIMEOUT}'"
    CENTAUR_LOG_REQUEST_FAILURES="-Dcentaur.log-request-failures=true"
    CENTAUR_CONF="${RUN_SPECIFIED_TEST_DIR_CMD} ${CENTAUR_LOG_REQUEST_FAILURES} ${CENTAUR_CROMWELL_MODE} ${CENTAUR_CROMWELL_JAR} ${CENTAUR_CROMWELL_CONF} ${CENTAUR_CROMWELL_LOG} ${CENTAUR_CROMWELL_RESTART} ${CENTAUR_SEND_RECEIVE_TIMEOUT}"
fi


TEST_DESCRIPTION="Running Centaur with sbt test"
TEST_COMMAND="java ${CENTAUR_CONF} -cp $CP org.scalatest.tools.Runner -R centaur/target/scala-2.12/it-classes -oD -u target/test-reports -PS${TEST_THREAD_COUNT}"

if [[ -n ${EXCLUDE_TAG[*]} ]]; then
    TEST_DESCRIPTION=${TEST_DESCRIPTION}" excluding ${EXCLUDE_TAG[*]} tests"
    EXCLUDE=""
    for val in "${EXCLUDE_TAG[@]}"; do
        EXCLUDE=" -l $val"${EXCLUDE}
    done
    TEST_COMMAND="${TEST_COMMAND}${EXCLUDE}"
fi

if [[ -n ${INCLUDE_TAG[*]} ]]; then
    TEST_DESCRIPTION=${TEST_DESCRIPTION}" including ${INCLUDE_TAG[*]} tests"
    INCLUDE=""
    for val in "${INCLUDE_TAG[@]}"; do
        INCLUDE=" -n $val"${INCLUDE}
    done
    TEST_COMMAND="${TEST_COMMAND}${INCLUDE}"
fi

if [[ -n "${SUITE}" ]]; then
    TEST_DESCRIPTION=${TEST_DESCRIPTION}" with suite ${SUITE}"
    TEST_COMMAND="${TEST_COMMAND} -s ${SUITE}"
fi

echo "${TEST_DESCRIPTION}"
echo "TEST_COMMAND is ${TEST_COMMAND}"

eval "${TEST_COMMAND} >> ${CENTAUR_LOG} 2>&1"

if [ $? -eq 0 ]; then
    EXIT_CODE=0
    TEST_STATUS="succeeded"
fi

echo "SBT test $TEST_STATUS, please see logs for more information"
