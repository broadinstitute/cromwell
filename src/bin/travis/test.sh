#!/usr/bin/env bash

set -e

SCRIPT_DIR=src/bin/travis
# $TRAVIS_EVENT_TYPE will be cron is this build was initiated by a cron job
if [ "$TRAVIS_EVENT_TYPE" == "cron" ]; then
    # BUILD_TYPE is coming in from the Travis build matrix
    if [ "$BUILD_TYPE" = "centaurJes" ]; then
        "${SCRIPT_DIR}"/testCentaurJes.sh
    elif [ "$BUILD_TYPE" = "centaurTes" ]; then
        "${SCRIPT_DIR}"/testCentaurTes.sh
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
else
    if [ "$BUILD_TYPE" = "centaurJes" ]; then
        "${SCRIPT_DIR}"/testCentaurJes.sh src/main/resources/integrationTestCases
    else
        exit 0
    fi
fi
