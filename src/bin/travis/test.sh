#!/usr/bin/env bash

set -e

SCRIPT_DIR=src/bin/travis

# BUILD_TYPE is coming in from the Travis build matrix
if [ "$BUILD_TYPE" = "centaurJes" ]; then
    "${SCRIPT_DIR}"/testCentaurJes.sh
elif [ "$BUILD_TYPE" = "centaurLocal" ]; then
    "${SCRIPT_DIR}"/testCentaurLocal.sh
elif [ "$BUILD_TYPE" = "sbt" ]; then
    "${SCRIPT_DIR}"/testSbt.sh
elif [ "$BUILD_TYPE" = "checkPublish" ]; then
    "${SCRIPT_DIR}"/testCheckPublish.sh
else
    echo "Unknown BUILD_TYPE: '$BUILD_TYPE'"
    exit 1
fi
