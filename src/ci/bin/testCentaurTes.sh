#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
export CROMWELL_BUILD_OPTIONAL_SECURE=true
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

cromwell::build::assemble_jars

FUNNEL_PATH="${CROMWELL_BUILD_ROOT_DIRECTORY}/funnel"
FUNNEL_CONF="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/funnel.conf"

# Increase max open files to the maximum allowed. Attempt to help on macos due to the default soft ulimit -n -S 256.
ulimit -n "$(ulimit -n -H)"
if [[ ! -f "${FUNNEL_PATH}" ]]; then
    FUNNEL_TAR_GZ="funnel-${CROMWELL_BUILD_OS}-amd64-0.5.0.tar.gz"
    curl "https://github.com/ohsu-comp-bio/funnel/releases/download/0.5.0/${FUNNEL_TAR_GZ}" -o "${FUNNEL_TAR_GZ}" -L
    tar xzf "${FUNNEL_TAR_GZ}"
fi

shutdown_funnel() {
    if [[ -n "${FUNNEL_PID+set}" ]]; then
        cromwell::build::kill_tree "${FUNNEL_PID}"
    fi
}

cromwell::build::add_exit_function shutdown_funnel

mkdir -p logs
nohup "${FUNNEL_PATH}" server run --config "${FUNNEL_CONF}" &> logs/funnel.log &

FUNNEL_PID=$!

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
