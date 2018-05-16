#!/usr/bin/env bash

# `sbt assembly` must have already been run.
BUILD_ROOT="$( dirname "${BASH_SOURCE[0]}" )/../../.."
CENTAUR_CWL_JAR="${CENTAUR_CWL_JAR:-"$( find "${BUILD_ROOT}/centaurCwlRunner/target/scala-2.12" -name 'centaur-cwl-runner-*.jar' | head -n 1 )"}"
CENTAUR_CWL_SKIP_FILE="${BUILD_ROOT}/centaurCwlRunner/src/main/resources/skipped_tests.csv"

java ${CENTAUR_CWL_JAVA_ARGS-"-Xmx1g"} -jar "${CENTAUR_CWL_JAR}" --skip-file "${CENTAUR_CWL_SKIP_FILE}" "$@"
