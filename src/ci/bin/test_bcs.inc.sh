#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail
# import in shellcheck / CI / IntelliJ compatible ways
# shellcheck source=/dev/null
source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of common BCS functions for use in other scripts.
#
# Functions:
#
#   - cromwell::build::bcs::*
#     Functions for use in other BCS scripts
#
#   - cromwell::private::bcs::*
#     Functions for use only within this file by cromwell::build::bcs::* functions
#

cromwell::private::bcs::bcs_install() {
    cromwell::build::pip_install batchcompute-cli==1.7.1 --upgrade
}

cromwell::private::bcs::bcs_run() {
    local bcs_command
    local bcs_command_result

    bcs_command="${1:?bcs_run called without a command}"; shift
    bcs_command_result="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${bcs_command}_result.txt"

    # Login failures print out the access key so send output to a file.
    bcs "${bcs_command}" "$@" </dev/null 2>&1 | tee "${bcs_command_result}"

    # bcs commands always exit with zero. make sure the result text does not contain "error".
    if grep -q -i error "${bcs_command_result}"; then
        echo "bcs ${bcs_command} failed" >&2
        grep -i error "${bcs_command_result}" >&2
        return 1
    else
        return 0
    fi
}

cromwell::private::bcs::try_bcs_login() {
    local bcs_login_include
    bcs_login_include="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/bcs_login.inc.sh"
    if [[ -f "${bcs_login_include}" ]]; then
        # shellcheck source=/dev/null
        source "${bcs_login_include}"
    fi
}

cromwell::private::bcs::try_bcs_create_cluster() {
    local cluster_name

    cluster_name=$(echo "cromwell_build_${CROMWELL_BUILD_PROVIDER}_${CROMWELL_BUILD_NUMBER}" | tr -c a-zA-Z0-9_- _)

    echo "Creating BCS cluster name: '${cluster_name}'"

    CROMWELL_BUILD_BCS_CLUSTER_ID=$( \
        cromwell::private::bcs::bcs_run \
            create_cluster \
            "${cluster_name}" \
            --image img-ubuntu \
            --type ecs.sn1ne.large \
            --nodes 8 \
            --vpc_cidr_block 192.168.1.0/24 \
            --no_cache_support \
            | grep 'Cluster created:' \
            | awk '{print $NF}' \
    )

    echo "Created BCS cluster id: '${CROMWELL_BUILD_BCS_CLUSTER_ID}'"

    export CROMWELL_BUILD_BCS_CLUSTER_ID
}

cromwell::private::bcs::try_bcs_delete_cluster() {
    cromwell::private::bcs::bcs_run delete_cluster --yes "${CROMWELL_BUILD_BCS_CLUSTER_ID}"
}

cromwell::private::bcs::bcs_login() {
    cromwell::build::exec_retry_function cromwell::private::bcs::try_bcs_login
}

cromwell::private::bcs::bcs_config() {
    cromwell::private::bcs::bcs_run config --god true
}

cromwell::private::bcs::bcs_create_cluster() {
    cromwell::build::exec_retry_function cromwell::private::bcs::try_bcs_create_cluster
    cromwell::build::add_exit_function cromwell::private::bcs::bcs_delete_cluster
}

cromwell::private::bcs::bcs_delete_cluster() {
    if [[ -n "${CROMWELL_BUILD_BCS_CLUSTER_ID}" ]]; then
        cromwell::build::exec_retry_function cromwell::private::bcs::try_bcs_delete_cluster || true
    fi
}

cromwell::private::bcs::bcs_delete_old_resources() {
    # Clean up assuming that all BCS jobs and clusters that are older than 3 hours are orphans to be deleted. jq 1.5
    # date functions all use UTC. https://stedolan.github.io/jq/manual/v1.5/#Dates Set the timezone environment variable
    # to UTC before running the command just in case this script/bcs are run on a different zone outside of UTC.
    echo "Please wait, removing old jobs…"

    TZ=utc \
        bcs job --all --show_json \
        | jq \
            -L "${CROMWELL_BUILD_RESOURCES_DIRECTORY}" \
            --raw-output 'include "bcs"; printIdsMoreThanSecondsOld(.; 3 * 60 * 60)' \
        | xargs -n 1 -I '{}' bash -c 'bcs delete_job --yes {} || true'

    echo "Please wait, removing old clusters…"
    TZ=utc \
        bcs cluster --show_json \
        | jq \
            -L "${CROMWELL_BUILD_RESOURCES_DIRECTORY}" \
            --raw-output 'include "bcs"; printIdsMoreThanSecondsOld(.; 3 * 60 * 60)' \
        | xargs -n 1 -I '{}' bash -c 'bcs delete_cluster --yes {} || true'
}

cromwell::build::bcs::setup_bcs_environment() {
    cromwell::private::bcs::bcs_install
    cromwell::build::exec_silent_function cromwell::private::bcs::bcs_login
    cromwell::private::bcs::bcs_config
    cromwell::private::bcs::bcs_delete_old_resources

    # Create the BCS cluster before sbt assembly as cluster creation takes a few minutes
    cromwell::private::bcs::bcs_create_cluster
}
