#!/bin/bash

set -e
set -x

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# BUILD_TYPE is coming in from the Travis build matrix
if [ "$BUILD_TYPE" = "centaur" ]; then
    "${SCRIPT_DIR}"/travisCentaur.sh
elif [ "$BUILD_TYPE" = "test" ]; then
    "${SCRIPT_DIR}"/travisTest.sh
else
    echo "I should not be here"
    exit 1
fi

