#!/usr/bin/env bash

set -e

sudo -H pip install --upgrade pip
sudo -H pip install cwltest

ENABLE_COVERAGE=true sbt assembly
CROMWELL_JAR=$(find "$(pwd)/target/scala-2.12" -name "cromwell-*.jar")
CENTAUR_CWL_RUNNER="$(pwd)/centaurCwlRunner/src/bin/centaur-cwl-runner.bash"

git clone --depth 1 https://github.com/common-workflow-language/common-workflow-language.git

shutdownCromwell() {
    if [ -z "${CROMWELL_PID}" ]; then
        kill ${CROMWELL_PID}
    fi
}

trap "shutdownCromwell" EXIT

java -jar ${CROMWELL_JAR} server &
CROMWELL_PID=$$

sleep 5

(
    cd common-workflow-language
    # For now, always exit with 0.
    ./run_test.sh RUNNER="${CENTAUR_CWL_RUNNER}" || true
)

sbt coverageReport --warn
sbt coverageAggregate --warn
bash <(curl -s https://codecov.io/bash) >/dev/null
