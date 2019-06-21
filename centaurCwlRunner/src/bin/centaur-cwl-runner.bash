#!/usr/bin/env bash

# `sbt assembly` must have already been run.
build_root="$( dirname "${BASH_SOURCE[0]}" )/../../.."
centaur_cwl_jar="${CENTAUR_CWL_JAR:-"$( find "${build_root}/centaurCwlRunner/target/scala-2.12" -name 'centaur-cwl-runner-*.jar' | head -n 1 )"}"
centaur_cwl_skip_file="${build_root}/centaurCwlRunner/src/main/resources/skipped_tests.csv"

centaur_cwl_java_args=("-Xmx1g")
if [[ -n "${CENTAUR_CWL_JAVA_ARGS-}" ]]; then
    # Allow splitting on space to simulate an exported array
    # https://stackoverflow.com/questions/5564418/exporting-an-array-in-bash-script#answer-5564589
    # shellcheck disable=SC2206
    centaur_cwl_java_args+=(${CENTAUR_CWL_JAVA_ARGS})
fi

# Handle empty arrays in older versions of bash
# https://stackoverflow.com/questions/7577052/bash-empty-array-expansion-with-set-u#answer-7577209
java \
    ${centaur_cwl_java_args[@]+"${centaur_cwl_java_args[@]}"} \
    -jar "${centaur_cwl_jar}" \
    --skip-file "${centaur_cwl_skip_file}" \
    "$@"
