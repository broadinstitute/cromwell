#!/usr/bin/env bash

set -e

SCRIPT_DIR=src/bin/travis

# Delete ~/.sbt/boot to fix consistent, almost immediate failures on sub-builds (usually TES but sometimes others).
# Even purging Travis caches didn't always fix the problem. Fortunately stackoverflow knew what to do:
# https://stackoverflow.com/questions/24539576/sbt-scala-2-10-4-missing-scala-tools-nsc-global
rm -rf ~/.sbt/boot/

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
        centaurBcs)
            "${SCRIPT_DIR}"/testCentaurBcs.sh
            ;;
        centaurCwlConformanceLocal)
            "${SCRIPT_DIR}"/testCentaurCwlConformanceLocal.sh
            ;;
        centaurCwlConformancePAPI)
            "${SCRIPT_DIR}"/testCentaurCwlConformancePAPI.sh
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
        centaurTes|centaurLocal|centaurBcs|centaurCwlConformanceLocal|centaurCwlConformancePAPI|sbt|checkPublish)
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
        exit 1
        ;;
    esac
