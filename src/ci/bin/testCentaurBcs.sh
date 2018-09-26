#!/usr/bin/env bash

# Disabled until:
# - https://github.com/broadinstitute/cromwell/issues/3555
# - https://github.com/broadinstitute/cromwell/issues/3554
exit 0

set -e
export CROMWELL_BUILD_REQUIRES_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

cromwell::build::setup_secure_resources

# https://github.com/broadinstitute/cromwell/issues/3522
# https://github.com/broadinstitute/cromwell/issues/3523
# https://github.com/broadinstitute/cromwell/issues/3524
exclude_known_bugs=( \
    -e bad_file_string \
    -e bad_output_task \
    -e tmp_dir \
)

# https://github.com/broadinstitute/cromwell/issues/3518
exclude_docker_tests=( \
    -e curl \
    -e docker_hash_dockerhub \
    -e docker_hash_gcr \
    -e docker_hash_quay \
    -e dont_cache_to_failed_jobs \
    -e hello \
    -e hello_yaml \
    -e inline_file \
    -e inline_file_custom_entryname \
    -e iwdr_input_string \
    -e iwdr_input_string_function \
    -e non_root_default_user \
    -e three_step_cwl \
)

# https://github.com/broadinstitute/cromwell/issues/3519
exclude_glob_tests=( \
    -e cwl_glob_sort \
    -e cwl_interpolated_strings \
    -e dontglobinputs \
    -e globbingbehavior \
    -e globbingindex \
    -e globbingscatter \
    -e lots_of_inputs \
    -e space \
    -e wdl_empty_glob \
)

centaur/test_cromwell.sh \
    -n "${CROMWELL_BUILD_CENTAUR_CONFIG}" \
    -l "${CROMWELL_BUILD_LOG_DIRECTORY}" \
    -g \
    -p 100 \
    -t 1m \
    -e localdockertest \
    "${exclude_known_bugs[@]}" \
    "${exclude_docker_tests[@]}" \
    "${exclude_glob_tests[@]}"

cromwell::build::generate_code_coverage
