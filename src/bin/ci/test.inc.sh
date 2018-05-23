#!/usr/bin/env bash

# A set of common functions for use in other scripts.
#
# Functions:
#
#   - cromwell::build::*
#     Functions for use in other scripts
#
#   - cromwell::private::*
#     Functions for use only within this file by cromwell::build::* functions
#
# Special Variables
#
#   - CROMWELL_BUILD_*
#     Variables for use in other scripts.
#
#   - CROMWELL_SECURE_*
#     Should not be used/printed/echoed!
#
#   - crmdbg
#     Quick debug scripts. Example: `crmdbg=y src/bin/ci/testCentaulLocal.sh`
#
#   - crmcrn
#     Simulate a cron build. Example: `crmcrn=y src/bin/ci/testCentaurPapiV2.sh`
#
cromwell::private::check_debug() {
    # shellcheck disable=SC2154
    if [ -n "${crmdbg:+set}" ]; then
        set -x
    fi

    # shellcheck disable=SC2154
    if [ -n "${crmcrn:+set}" ]; then
        CROMWELL_BUILD_IS_CRON=true
    fi
}

# Exports environment variables used for scripts.
cromwell::private::create_build_variables() {
    if [ "${TRAVIS}" = "true" ]; then
        CROMWELL_BUILD_IS_CI=true
    else
        CROMWELL_BUILD_IS_CI=false
    fi

    if [ -n "${CROMWELL_BUILD_IS_CRON}" ]; then
        # ignore
        :
    elif [ "${TRAVIS_EVENT_TYPE}" = "cron" ]; then
        CROMWELL_BUILD_IS_CRON=true
    else
        CROMWELL_BUILD_IS_CRON=false
    fi

    if [ "${TRAVIS_SECURE_ENV_VARS}" = "true" ]; then
        CROMWELL_BUILD_IS_SECURE=true
    else
        CROMWELL_BUILD_IS_SECURE=false
    fi

    CROMWELL_BUILD_TYPE="${BUILD_TYPE}"
    CROMWELL_BUILD_BRANCH="${TRAVIS_BRANCH}"
    CROMWELL_BUILD_EVENT="${TRAVIS_EVENT_TYPE}"
    CROMWELL_BUILD_TAG="${TRAVIS_TAG}"

    # simplified from https://stackoverflow.com/a/18434831/3320205
    CROMWELL_BUILD_OS_DARWIN="darwin";
    CROMWELL_BUILD_OS_LINUX="linux";
    case "${OSTYPE}" in
        darwin*)  CROMWELL_BUILD_OS="${CROMWELL_BUILD_OS_DARWIN}" ;;
        linux*)   CROMWELL_BUILD_OS="${CROMWELL_BUILD_OS_LINUX}" ;;
        *)        CROMWELL_BUILD_OS="unknown_os" ;;
    esac

    CROMWELL_BUILD_HOME_DIRECTORY="${HOME}"
    CROMWELL_BUILD_ROOT_DIRECTORY="$(pwd)"

    CROMWELL_BUILD_SCRIPTS_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/src/bin/ci"
    CROMWELL_BUILD_SCRIPTS_RESOURCES="${CROMWELL_BUILD_SCRIPTS_DIRECTORY}/resources"

    CROMWELL_BUILD_GIT_USER_EMAIL="travis@travis-ci.org"
    CROMWELL_BUILD_GIT_USER_NAME="Travis CI"

    CROMWELL_BUILD_MYSQL_USERNAME="root"
    CROMWELL_BUILD_CENTAUR_RESOURCES="${CROMWELL_BUILD_ROOT_DIRECTORY}/centaur/src/main/resources"
    CROMWELL_BUILD_CENTAUR_INTEGRATION_TESTS="${CROMWELL_BUILD_CENTAUR_RESOURCES}/integrationTestCases"

    CROMWELL_BUILD_EXIT_FUNCTIONS="${CROMWELL_BUILD_ROOT_DIRECTORY}/cromwell_build_exit_functions.$$"

    export CROMWELL_BUILD_IS_CI
    export CROMWELL_BUILD_IS_CRON
    export CROMWELL_BUILD_IS_SECURE
    export CROMWELL_BUILD_TYPE
    export CROMWELL_BUILD_BRANCH
    export CROMWELL_BUILD_EVENT
    export CROMWELL_BUILD_TAG
    export CROMWELL_BUILD_OS
    export CROMWELL_BUILD_OS_DARWIN
    export CROMWELL_BUILD_OS_LINUX
    export CROMWELL_BUILD_HOME_DIRECTORY
    export CROMWELL_BUILD_ROOT_DIRECTORY
    export CROMWELL_BUILD_SCRIPTS_DIRECTORY
    export CROMWELL_BUILD_SCRIPTS_RESOURCES
    export CROMWELL_BUILD_GIT_USER_EMAIL
    export CROMWELL_BUILD_GIT_USER_NAME
    export CROMWELL_BUILD_MYSQL_USERNAME
    export CROMWELL_BUILD_CENTAUR_RESOURCES
    export CROMWELL_BUILD_CENTAUR_INTEGRATION_TESTS
    export CROMWELL_BUILD_EXIT_FUNCTIONS
}

cromwell::private::echo_build_variables() {
    echo "CROMWELL_BUILD_TYPE='${CROMWELL_BUILD_TYPE}'"
    echo "CROMWELL_BUILD_BRANCH='${CROMWELL_BUILD_BRANCH}'"
    echo "CROMWELL_BUILD_EVENT='${CROMWELL_BUILD_EVENT}'"
    echo "CROMWELL_BUILD_TAG='${CROMWELL_BUILD_TAG}'"
}

cromwell::private::verify_secure_variables() {
    if [ "${CROMWELL_BUILD_IS_SECURE}" = "false" ]; then
        echo "********************************************************"
        echo "********************************************************"
        echo "**                                                    **"
        echo "**  WARNING: Encrypted keys are unavailable. Exiting. **"
        echo "**                                                    **"
        echo "********************************************************"
        echo "********************************************************"
        exit 0
    fi
}

cromwell::private::verify_cron_build() {
    if [ "${CROMWELL_BUILD_IS_CRON}" = "true" ] && [ "${CROMWELL_BUILD_SUPPORTS_CRON}" != "true" ]; then
        echo "***************************"
        echo "***************************"
        echo "**                       **"
        echo "**  Skipping cron build. **"
        echo "**                       **"
        echo "***************************"
        echo "***************************"
        exit 0
    fi
}

cromwell::private::setup_secure_variables() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::setup_secure_variables
    else
        CROMWELL_SECURE_DOCKER_USERNAME="${DOCKER_USERNAME}"
        CROMWELL_SECURE_DOCKER_PASSWORD="${DOCKER_PASSWORD}"
        CROMWELL_SECURE_VAULT_TOKEN="${JES_TOKEN}"
    fi
}

cromwell::private::export_conformance_variables() {
    CROMWELL_BUILD_CWL_TEST_RUNNER="${CROMWELL_BUILD_ROOT_DIRECTORY}/centaurCwlRunner/src/bin/centaur-cwl-runner.bash"
    CROMWELL_BUILD_CWL_TEST_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/common-workflow-language"
    CROMWELL_BUILD_CWL_TEST_RESOURCES="${CROMWELL_BUILD_CWL_TEST_DIRECTORY}/v1.0/v1.0"
    CROMWELL_BUILD_CWL_TEST_WDL="${CROMWELL_BUILD_SCRIPTS_RESOURCES}/cwl_conformance_test.wdl"
    CROMWELL_BUILD_CWL_TEST_INPUTS="${CROMWELL_BUILD_SCRIPTS_RESOURCES}/cwl_conformance_test.inputs.json"
    CROMWELL_BUILD_CWL_TEST_OUTPUT="${CROMWELL_BUILD_ROOT_DIRECTORY}/cwl_conformance_test.out.txt"
    CROMWELL_BUILD_CWL_TEST_COMMIT=e7ae597e79d426feca3d4faa474806070c3fa7d5 # use known git hash to avoid tests changes
    CROMWELL_BUILD_CWL_TEST_PARALLELISM=10 # Set too high will cause false negatives due to cromwell server timeouts.

    export CROMWELL_BUILD_CWL_TEST_RUNNER
    export CROMWELL_BUILD_CWL_TEST_DIRECTORY
    export CROMWELL_BUILD_CWL_TEST_RESOURCES
    export CROMWELL_BUILD_CWL_TEST_WDL
    export CROMWELL_BUILD_CWL_TEST_INPUTS
    export CROMWELL_BUILD_CWL_TEST_OUTPUT
    export CROMWELL_BUILD_CWL_TEST_COMMIT
    export CROMWELL_BUILD_CWL_TEST_PARALLELISM
}

cromwell::private::exec_test_script() {
    local upper_build_type
    upper_build_type="$(tr '[:lower:]' '[:upper:]' <<< "${CROMWELL_BUILD_TYPE:0:1}")${CROMWELL_BUILD_TYPE:1}"
    exec "${CROMWELL_BUILD_SCRIPTS_DIRECTORY}/test${upper_build_type}.sh"
}

cromwell::private::delete_boto_config() {
    # https://github.com/travis-ci/travis-ci/issues/7940#issuecomment-310759657
    sudo rm -f /etc/boto.cfg
    export BOTO_CONFIG=/dev/null
}

cromwell::private::delete_sbt_boot() {
    # Delete ~/.sbt/boot to fix consistent, almost immediate failures on sub-builds (usually TES but sometimes others).
    # Even purging Travis caches didn't always fix the problem. Fortunately stackoverflow knew what to do:
    # https://stackoverflow.com/questions/24539576/sbt-scala-2-10-4-missing-scala-tools-nsc-global
    rm -rf ~/.sbt/boot/
}

cromwell::private::upgrade_pip() {
    sudo -H pip install --upgrade pip
}

cromwell::private::create_mysql_cromwell_test() {
    mysql -u "${CROMWELL_BUILD_MYSQL_USERNAME}" -e "SET GLOBAL sql_mode = 'STRICT_ALL_TABLES';"
    mysql -u "${CROMWELL_BUILD_MYSQL_USERNAME}" -e "CREATE DATABASE IF NOT EXISTS cromwell_test;"
}

cromwell::private::pull_common_docker_images() {
    # All tests use ubuntu:latest - make sure it's there before starting the tests
    # because pulling the image during some of the tests would cause them to fail
    # (specifically output_redirection which expects a specific value in stderr)
    docker pull ubuntu
}

cromwell::private::install_cwltool() {
    # TODO: No clue why these are needed for cwltool. If you know please update this comment.
    sudo apt-get install procps || true
    sudo -H pip install 'requests[security]'
    sudo -H pip install --ignore-installed cwltool
}

cromwell::private::install_cwltest() {
    sudo -H pip install cwltest
}

cromwell::private::checkout_pinned_cwl() {
    if [ ! -d "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}" ]; then
        git clone \
            https://github.com/common-workflow-language/common-workflow-language.git \
            "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}"
        (
            cd "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}" || exit 1
            git checkout "${CROMWELL_BUILD_CWL_TEST_COMMIT}"
        )
    fi
}

cromwell::private::docker_login() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::docker_login
    else
        local dockerhub_auth_include
        dockerhub_auth_include="${CROMWELL_BUILD_SCRIPTS_RESOURCES}/dockerhub_auth.sh"
        if [ -f "${dockerhub_auth_include}" ]; then
            # shellcheck source=/dev/null
            source "${dockerhub_auth_include}"
        fi
        docker login \
            --username "${CROMWELL_SECURE_DOCKER_USERNAME}" \
            --password "${CROMWELL_SECURE_DOCKER_PASSWORD}"
    fi
}

cromwell::private::vault_login() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::vault_login
    else
        # Login to vault to access secrets
        docker run --rm \
            -v "${CROMWELL_BUILD_HOME_DIRECTORY}:/root:rw" \
            broadinstitute/dsde-toolbox \
            vault auth "${CROMWELL_SECURE_VAULT_TOKEN}" < /dev/null > /dev/null && echo vault auth success
    fi
}

cromwell::private::render_secure_resources() {
    docker run --rm \
        -v "${CROMWELL_BUILD_HOME_DIRECTORY}:/root:rw" \
        -v "${CROMWELL_BUILD_SCRIPTS_RESOURCES}:/resources" \
        -e ENVIRONMENT=not_used \
        -e INPUT_PATH=/resources \
        -e OUT_PATH=/resources \
        broadinstitute/dsde-toolbox render-templates.sh \
    || {
        echo "**************************************************************"
        echo "**************************************************************"
        echo "**                                                          **"
        echo "**        WARNING: Unable to render vault resources.        **"
        echo "**  '*.ctmpl' files should be copied and updated manually.  **"
        echo "**                                                          **"
        echo "**************************************************************"
        echo "**************************************************************"
    }
}

cromwell::private::setup_secure_resources() {
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
            cromwell::private::setup_secure_variables
            # Ignore premature login errors that occur once dsde-toolbox is public and u/p are only in vault
            cromwell::private::docker_login || true
            cromwell::private::vault_login
            cromwell::private::render_secure_resources
            cromwell::private::docker_login
    else
            cromwell::private::render_secure_resources
    fi
}

cromwell::private::generate_code_coverage() {
    sbt coverageReport -warn
    sbt coverageAggregate -warn
    bash <(curl -s https://codecov.io/bash) > /dev/null || true
}

cromwell::private::publish_artifacts_only() {
    sbt "$@" +publish
}

cromwell::private::publish_artifacts_and_docker() {
    sbt "$@" +publish dockerBuildAndPush
}

# Some CI environments want to know when new docker images are published. They do not currently poll dockerhub but do
# poll github. To help those environments, signal that a new set of docker images has been published to dockerhub by
# updating a well known branch in github.
cromwell::private::push_publish_complete() {
    local github_private_deploy_key="${CROMWELL_BUILD_SCRIPTS_RESOURCES}/github_private_deploy_key"
    local git_repo="git@github.com:broadinstitute/cromwell.git"
    local git_branch="${CROMWELL_BUILD_BRANCH}_publish_complete"
    local git_remote="publish_complete"
    local git_message="publish complete [skip ci]"

    # Loosely adapted from https://github.com/broadinstitute/workbench-libs/blob/435a932/scripts/version_update.sh
    mkdir publish_complete
    (
        cd publish_complete || exit 1

        git init
        git config core.sshCommand "ssh -i $github_private_deploy_key -F /dev/null"
        git config user.email "${CROMWELL_BUILD_GIT_USER_EMAIL}"
        git config user.name "${CROMWELL_BUILD_GIT_USER_NAME}"

        git remote add "$git_remote" "$git_repo"
        git checkout -b "$git_branch"
        git commit --allow-empty -m "$git_message"
        git push -f "$git_remote" "$git_branch"
    )
}

cromwell::private::start_build_heartbeat() {
    # Sleep one minute between printouts, but don't zombie for more than three hours
    for ((i=0; i < 180; i++)); do
        sleep 60
        printf "…"
    done &
    CROMWELL_BUILD_HEARTBEAT_PID=$!
}

cromwell::private::start_cromwell_log_tail() {
    while [ ! -f logs/cromwell.log ]; do
        sleep 2
    done && tail -n 0 -f logs/cromwell.log 2> /dev/null &
    CROMWELL_BUILD_CROMWELL_LOG_TAIL_PID=$!
}

cromwell::private::start_centaur_log_tail() {
    while [ ! -f logs/centaur.log ]; do
        sleep 2
    done && tail -n 0 -f logs/centaur.log 2> /dev/null &
    CROMWELL_BUILD_CENTAUR_LOG_TAIL_PID=$!
}

cromwell::private::cat_centaur_log() {
    echo "CENTAUR LOG"
    cat logs/centaur.log
}

cromwell::private::cat_conformance_log() {
    echo "CONFORMANCE LOG"
    cat "${CROMWELL_BUILD_CWL_TEST_OUTPUT}"
}

cromwell::private::kill_build_heartbeat() {
    if [ -n "${CROMWELL_BUILD_HEARTBEAT_PID:+set}" ]; then
        cromwell::private::kill_tree "${CROMWELL_BUILD_HEARTBEAT_PID}"
    fi
}

cromwell::private::kill_cromwell_log_tail() {
    if [ -n "${CROMWELL_BUILD_CROMWELL_LOG_TAIL_PID:+set}" ]; then
        cromwell::private::kill_tree "${CROMWELL_BUILD_CROMWELL_LOG_TAIL_PID}"
    fi
}

cromwell::private::kill_centaur_log_tail() {
    if [ -n "${CROMWELL_BUILD_CENTAUR_LOG_TAIL_PID:+set}" ]; then
        cromwell::private::kill_tree ${CROMWELL_BUILD_CENTAUR_LOG_TAIL_PID}
    fi
}

cromwell::private::run_exit_functions() {
    if [ -f "${CROMWELL_BUILD_EXIT_FUNCTIONS}" ]; then
        local exit_function
        while read exit_function; do
          ${exit_function} || true
        done < "${CROMWELL_BUILD_EXIT_FUNCTIONS}"
        rm "${CROMWELL_BUILD_EXIT_FUNCTIONS}" || true
    fi
}

# Adds the function to the list of functions to run on exit.
# Requires one positional parameter, the function to run.
cromwell::private::add_exit_function() {
    local exit_function
    exit_function="${1:?add_exit_function called without a function}"
    echo "${exit_function}" >> "${CROMWELL_BUILD_EXIT_FUNCTIONS}"
    trap cromwell::private::run_exit_functions TERM EXIT
}

cromwell::private::exec_silent_function() {
    local silent_function
    local xtrace_restore_function
    silent_function="${1:?exec_silent_function called without a function}"
    xtrace_restore_function="$(shopt -po xtrace)"
    shopt -uo xtrace
    ${silent_function}
    ${xtrace_restore_function}
}

cromwell::private::is_xtrace_enabled() {
    # Return 0 if xtrace is disabled (set +x), 1 if xtrace is enabled (set -x).
    shopt -qo xtrace
}

cromwell::private::kill_tree() {
  local pid
  local cpid
  pid="${1:?kill_tree called without a pid}"
  for cpid in $(pgrep -P "${pid}"); do
    cromwell::private::kill_tree "${cpid}"
  done
  kill "${pid}" 2> /dev/null
}

cromwell::build::exec_test_script() {
    cromwell::private::create_build_variables
    cromwell::private::echo_build_variables
    cromwell::private::exec_test_script
}

cromwell::build::setup_common_environment() {
    cromwell::private::check_debug
    cromwell::private::create_build_variables
    cromwell::private::verify_cron_build
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
        cromwell::private::delete_boto_config
        cromwell::private::delete_sbt_boot
        cromwell::private::upgrade_pip
        cromwell::private::create_mysql_cromwell_test
        cromwell::private::pull_common_docker_images
        cromwell::private::install_cwltool
    fi
}

cromwell::build::setup_secure_resources() {
    cromwell::private::setup_secure_resources
}

cromwell::build::setup_centaur_environment() {
    cromwell::private::start_build_heartbeat
    cromwell::private::start_cromwell_log_tail
    cromwell::private::start_centaur_log_tail
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
        cromwell::private::add_exit_function cromwell::private::cat_centaur_log
    fi
    cromwell::private::add_exit_function cromwell::private::kill_build_heartbeat
    cromwell::private::add_exit_function cromwell::private::kill_cromwell_log_tail
    cromwell::private::add_exit_function cromwell::private::kill_centaur_log_tail
}

cromwell::build::setup_conformance_environment() {
    cromwell::private::export_conformance_variables
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
        cromwell::private::install_cwltest
    fi
    cromwell::private::checkout_pinned_cwl
    cromwell::private::start_build_heartbeat
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
        cromwell::private::add_exit_function cromwell::private::cat_conformance_log
    fi
    cromwell::private::add_exit_function cromwell::private::kill_build_heartbeat
}

cromwell::build::assemble_jars() {
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
        CROMWELL_SBT_ASSEMBLY_LOG_LEVEL=error CROMWELL_SBT_COVERAGE=true sbt assembly -error
    fi
    CROMWELL_BUILD_JAR="$( \
        find "${CROMWELL_BUILD_ROOT_DIRECTORY}/server/target/scala-2.12" -name "cromwell-*.jar" \
        | head -n 1 \
        )"
    [ -s "${CROMWELL_BUILD_JAR}" ] || exit 1
    export CROMWELL_BUILD_JAR
}

cromwell::build::generate_code_coverage() {
    if [ "${CROMWELL_BUILD_IS_CRON}" != "true" ]; then
        cromwell::private::generate_code_coverage
    fi
}

cromwell::build::publish_artifacts() {
    if [ "${CROMWELL_BUILD_TYPE}" == "sbt" ] && [ "${CROMWELL_BUILD_EVENT}" == "push" ]; then

        if [ "${CROMWELL_BUILD_BRANCH}" = "develop" ]; then
            # Publish images for both the "cromwell develop branch" and the "cromwell dev environment".
            cromwell::private::setup_secure_resources
            CROMWELL_SBT_DOCKER_TAGS=develop,dev \
                cromwell::private::publish_artifacts_and_docker \
                -Dproject.isSnapshot=true
            cromwell::private::push_publish_complete
        fi

        if [[ "${CROMWELL_BUILD_BRANCH}" =~ ^[0-9\.]+_hotfix$ ]]; then
            # Docker tags float. "30" is the latest hotfix. Those dockers are published here on each hotfix commit.
            cromwell::private::setup_secure_resources
            cromwell::private::publish_artifacts_and_docker -Dproject.isSnapshot=false
        fi

        if [ -n "${CROMWELL_BUILD_TAG:+set}" ]; then
            # Artifact tags are static. Once "30" is set that is only "30" forever. Those artifacts are published here.
            cromwell::private::setup_secure_resources
            cromwell::private::publish_artifacts_only \
                -Dproject.version="${CROMWELL_BUILD_TAG}" \
                -Dproject.isSnapshot=false

        fi

    fi
}

cromwell::build::add_exit_function() {
    cromwell::private::add_exit_function "$1"
}

cromwell::build::kill_tree() {
    cromwell::private::kill_tree "$1"
}
