#!/bin/bash

shutdown() {
    cd "${INITIAL_DIR}"
    # This will take out the backgrounded Cromwell instance
    pkill -P $$
    exit "${EXIT_CODE}"
}

INITIAL_DIR=$(pwd)
EXIT_CODE=1

trap "shutdown" EXIT

set -e

# Used by the travis builds to build & run centaur against a prebuilt Cromwell. There is almost no
# handholding here so running it outside of Travis is likely to not do what you want
PARALLELISM_FACTOR=5

echo "Launching Cromwell for Centaur testing"
java -Dconfig.file=./jes.conf -jar "${INITIAL_DIR}"/target/scala-2.11/cromwell-*.jar server > cromwell.log &

# Build and run centaur
echo "Cloning Centaur"
git clone https://github.com/broadinstitute/centaur.git

cd centaur

echo "Running Centaur with ${PARALLELISM_FACTOR}-way parallelism"
./run_tests_parallel.sh "${PARALLELISM_FACTOR}"

if [ $? -eq 0 ]; then
  EXIT_CODE=0
fi
