#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

startup_funnel() {
    local funnel_path
    local funnel_conf
    local funnel_tar_gz

    funnel_path="${CROMWELL_BUILD_ROOT_DIRECTORY}/funnel"
    funnel_conf="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/funnel.conf"

    # Increase max open files to the maximum allowed. Attempt to help on macos due to the default soft ulimit -n -S 256.
    ulimit -n "$(ulimit -n -H)"
    if [[ ! -f "${funnel_path}" ]]; then
        funnel_tar_gz="funnel-${CROMWELL_BUILD_OS}-amd64-0.5.0.tar.gz"
        curl \
            --location \
            --output "${funnel_tar_gz}" \
            "https://github.com/ohsu-comp-bio/funnel/releases/download/0.5.0/${funnel_tar_gz}"
        tar xzf "${funnel_tar_gz}"
    fi

    mkdir -p logs
    nohup "${funnel_path}" server run --config "${funnel_conf}" &> logs/funnel.log &

    FUNNEL_PID=$!
    export FUNNEL_PID

    cromwell::build::add_exit_function shutdown_funnel
}

shutdown_funnel() {
    if [[ -n "${FUNNEL_PID+set}" ]]; then
        cromwell::build::kill_tree "${FUNNEL_PID}"
    fi
}

startup_funnel

# The following tests are skipped:
#
# call_cache_capoeira_local: fails on task 'read_files_without_docker' since the 'docker' runtime key is required for this backend
# draft3_call_cache_capoeira_local: same as above
# lots_of_inputs:            Funnel mounts in each input separately, this task surpasses the docker limit for volumes
# no_new_calls:              TES does not support checking job status after restart, and cannot tell if shouldSucceed is done or failed
# non_root_specified_user:   TES doesn't support switching users in the image
# write_lines_files:         all inputs are read-only in TES
# read_file_limits:          Fail only in Travis for unknown reason (Note that the draft 3 version does not fail)

cromwell::build::run_centaur \
    -e call_cache_capoeira_local \
    -e draft3_call_cache_capoeira_local \
    -e read_file_limits \
    -e lots_of_inputs \
    -e no_new_calls \
    -e non_root_default_user \
    -e non_root_specified_user \
    -e write_lines_files \
    -e draft3_read_write_functions_local \
    -e cwl_input_json \
    -e directory_type_local \

cromwell::build::generate_code_coverage
