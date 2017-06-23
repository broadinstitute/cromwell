#!/usr/bin/env bash

set -e

SCRIPT_DIR=src/bin/travis

if [ "$TRAVIS_EVENT_TYPE" = "push" ] || [ "$TRAVIS_EVENT_TYPE" = "pull_request" ] || [ "$TRAVIS_EVENT_TYPE" = "api" ]; then
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
fi

# $TRAVIS_EVENT_TYPE will be cron if this build was initiated by a cron job
if [ "$TRAVIS_EVENT_TYPE" = "cron" ]; then
    if [ "$BUILD_TYPE" = "centaurJes" ]; then
        "${SCRIPT_DIR}"/testCentaurJes.sh -i
    elif [ "$BUILD_TYPE" = "centaurTes" ] || [ "$BUILD_TYPE" = "centaurLocal" ] || [ "$BUILD_TYPE" = "sbt" ] || [ "$BUILD_TYPE" = "checkPublish" ]; then
        exit 0
    else
        echo "Unknown BUILD_TYPE: '$BUILD_TYPE'"
        exit 1
    fi
fi
