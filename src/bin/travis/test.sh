#!/usr/bin/env bash

set -e

SCRIPT_DIR=src/bin/travis

# $TRAVIS_EVENT_TYPE will be cron if this build was initiated by a cron job
case "$TRAVIS_EVENT_TYPE" in
    push|pull_request|api)
    # BUILD_TYPE is coming in from the Travis build matrix
        case "$BUILD_TYPE" in
        centaurJes)
            "${SCRIPT_DIR}"/testCentaurJes.sh
            ;;
        centaurTes)
            "${SCRIPT_DIR}"/testCentaurTes.sh
            ;;
        centaurLocal)
            "${SCRIPT_DIR}"/testCentaurLocal.sh
            ;;
        sbt)
            "${SCRIPT_DIR}"/testSbt.sh
            ;;
        checkPublish)
            "${SCRIPT_DIR}"/testCheckPublish.sh
            ;;
        *)
            echo "Unknown BUILD_TYPE: '$BUILD_TYPE'"
            exit 1
            ;;
        esac
        ;;
    cron)
        case "$BUILD_TYPE" in
        centaurJes)
            "${SCRIPT_DIR}"/testCentaurJes.sh -i
            ;;
        centaurTes|centaurLocal|sbt|checkPublish)
            exit 0
            ;;
        *)
            echo "Unknown BUILD_TYPE: '$BUILD_TYPE'"
            exit 1
            ;;
        esac
        ;;
    *)
        echo "Unknown TRAVIS_EVENT_TYPE: '$TRAVIS_EVENT_TYPE'"
    esac
