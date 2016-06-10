#!/bin/bash

shutdown() {
    cd "${INITIAL_DIR}"
    # This will take out the backgrounded Cromwell instance
    pkill -P $$
    exit 0
}

INITIAL_DIR=$(pwd)

trap "shutdown" EXIT

set -e

# Used by the travis builds to build & run centaur against a prebuilt Cromwell. There is almost no
# handholding here so running it outside of Travis is likely to not do what you want
PARALLELISM_FACTOR=3

java -jar ${INITIAL_DIR}/target/scala-2.11/cromwell-*.jar server > cromwell.log &

# Build and run centaur
echo "Cloning Centaur"
git clone https://github.com/broadinstitute/centaur.git

cd centaur

echo "Running Centaur with ${PARALLELISM_FACTOR}-way parallelism"
./run_tests_parallel.sh "${PARALLELISM_FACTOR}"

exit 0
