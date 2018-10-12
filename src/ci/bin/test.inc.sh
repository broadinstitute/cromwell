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
#   - crmdbg
#     Quick debug scripts. Example: `crmdbg=y src/ci/bin/testCentaulLocal.sh`
#
#   - crmcit
#     Simulate a centaur integration test build. Example: `crmcit=y src/ci/bin/testCentaurPapiV2.sh`
#
cromwell::private::check_debug() {
    # shellcheck disable=SC2154
    if [ -n "${crmdbg:+set}" ]; then
        set -x
    fi

    # shellcheck disable=SC2154
    if [ -n "${crmcit:+set}" ]; then
        CROMWELL_BUILD_CENTAUR_TYPE="integration"
    fi
}

# Exports environment variables used for scripts.
cromwell::private::create_build_variables() {
    CROMWELL_BUILD_PROVIDER_TRAVIS="travis"
    CROMWELL_BUILD_PROVIDER_JENKINS="jenkins"
    CROMWELL_BUILD_PROVIDER_UNKNOWN="unknown"

    if [ "${TRAVIS}" = "true" ]; then
        CROMWELL_BUILD_PROVIDER="${CROMWELL_BUILD_PROVIDER_TRAVIS}"
    elif [ "${JENKINS}" = "true" ]; then
        CROMWELL_BUILD_PROVIDER="${CROMWELL_BUILD_PROVIDER_JENKINS}"
    else
        CROMWELL_BUILD_PROVIDER="${CROMWELL_BUILD_PROVIDER_UNKNOWN}"
    fi

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
    CROMWELL_BUILD_LOG_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/target/ci/logs"
    CROMWELL_BUILD_CROMWELL_LOG="${CROMWELL_BUILD_LOG_DIRECTORY}/cromwell.log"

    CROMWELL_BUILD_SCRIPTS_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/src/ci/bin"
    CROMWELL_BUILD_RESOURCES_SOURCES="${CROMWELL_BUILD_ROOT_DIRECTORY}/src/ci/resources"
    CROMWELL_BUILD_RESOURCES_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/target/ci/resources"

    CROMWELL_BUILD_WAIT_FOR_IT_FILENAME="wait-for-it.sh"
    CROMWELL_BUILD_WAIT_FOR_IT_BRANCH="db049716e42767d39961e95dd9696103dca813f1"
    CROMWELL_BUILD_WAIT_FOR_IT_URL="https://raw.githubusercontent.com/vishnubob/wait-for-it/${CROMWELL_BUILD_WAIT_FOR_IT_BRANCH}/${CROMWELL_BUILD_WAIT_FOR_IT_FILENAME}"
    CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${CROMWELL_BUILD_WAIT_FOR_IT_FILENAME}"
    CROMWELL_BUILD_EXIT_FUNCTIONS="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_exit_functions.$$"

    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
            CROMWELL_BUILD_IS_CI=true
            CROMWELL_BUILD_IS_SECURE="${TRAVIS_SECURE_ENV_VARS}"
            CROMWELL_BUILD_TYPE="${BUILD_TYPE}"
            CROMWELL_BUILD_BRANCH="${TRAVIS_PULL_REQUEST_BRANCH:-${TRAVIS_BRANCH}}"
            CROMWELL_BUILD_EVENT="${TRAVIS_EVENT_TYPE}"
            CROMWELL_BUILD_TAG="${TRAVIS_TAG}"
            CROMWELL_BUILD_NUMBER="${TRAVIS_JOB_NUMBER}"
            CROMWELL_BUILD_URL="https://travis-ci.org/${TRAVIS_REPO_SLUG}/jobs/${TRAVIS_JOB_ID}"
            CROMWELL_BUILD_GIT_USER_EMAIL="travis@travis-ci.org"
            CROMWELL_BUILD_GIT_USER_NAME="Travis CI"
            CROMWELL_BUILD_VAULT_TOKEN="${CROMWELL_BUILD_HOME_DIRECTORY}/.vault-token"
            CROMWELL_BUILD_HEARTBEAT_MESSAGE="…"
            CROMWELL_BUILD_MYSQL_HOSTNAME="localhost"
            CROMWELL_BUILD_MYSQL_PORT="3306"
            CROMWELL_BUILD_MYSQL_USERNAME="travis"
            CROMWELL_BUILD_MYSQL_PASSWORD=""
            CROMWELL_BUILD_MYSQL_SCHEMA="cromwell_test"
            CROMWELL_BUILD_GENERATE_COVERAGE=true
            ;;
        "${CROMWELL_BUILD_PROVIDER_JENKINS}")
            # External variables must be passed through in the ENVIRONMENT of src/ci/docker-compose/docker-compose.yml
            CROMWELL_BUILD_IS_CI=true
            CROMWELL_BUILD_IS_SECURE=true
            CROMWELL_BUILD_TYPE="${JENKINS_BUILD_TYPE}"
            CROMWELL_BUILD_BRANCH="${GIT_BRANCH#origin/}"
            CROMWELL_BUILD_EVENT=""
            CROMWELL_BUILD_TAG=""
            CROMWELL_BUILD_NUMBER="${BUILD_NUMBER}"
            CROMWELL_BUILD_URL="${BUILD_URL}"
            CROMWELL_BUILD_GIT_USER_EMAIL="jenkins@jenkins.io"
            CROMWELL_BUILD_GIT_USER_NAME="Jenkins CI"
            CROMWELL_BUILD_VAULT_TOKEN="/dev/null"
            CROMWELL_BUILD_HEARTBEAT_MESSAGE="…\n"
            CROMWELL_BUILD_MYSQL_HOSTNAME="mysql-db"
            CROMWELL_BUILD_MYSQL_PORT="3306"
            CROMWELL_BUILD_MYSQL_USERNAME="root"
            CROMWELL_BUILD_MYSQL_PASSWORD=""
            CROMWELL_BUILD_MYSQL_SCHEMA="cromwell_test"
            CROMWELL_BUILD_GENERATE_COVERAGE=false
            ;;
        *)
            CROMWELL_BUILD_IS_CI=false
            CROMWELL_BUILD_IS_SECURE=true
            CROMWELL_BUILD_TYPE="unknown"
            CROMWELL_BUILD_BRANCH="unknown"
            CROMWELL_BUILD_EVENT="unknown"
            CROMWELL_BUILD_TAG=""
            CROMWELL_BUILD_NUMBER=""
            CROMWELL_BUILD_URL=""
            CROMWELL_BUILD_GIT_USER_EMAIL="unknown.git.user@example.org"
            CROMWELL_BUILD_GIT_USER_NAME="Unknown Git User"
            CROMWELL_BUILD_VAULT_TOKEN="${CROMWELL_BUILD_HOME_DIRECTORY}/.vault-token"
            CROMWELL_BUILD_HEARTBEAT_MESSAGE="…"
            CROMWELL_BUILD_MYSQL_HOSTNAME="${CROMWELL_BUILD_MYSQL_HOSTNAME-localhost}"
            CROMWELL_BUILD_MYSQL_PORT="${CROMWELL_BUILD_MYSQL_PORT-3306}"
            CROMWELL_BUILD_MYSQL_USERNAME="${CROMWELL_BUILD_MYSQL_USERNAME-root}"
            CROMWELL_BUILD_MYSQL_PASSWORD="${CROMWELL_BUILD_MYSQL_PASSWORD-}"
            CROMWELL_BUILD_MYSQL_SCHEMA="${CROMWELL_BUILD_MYSQL_SCHEMA-cromwell_test}"
            CROMWELL_BUILD_GENERATE_COVERAGE=true

            local bash_script
            for bash_script in "${BASH_SOURCE[@]}"; do
                if [[ "${bash_script}" != */test.inc.sh ]]; then
                    local build_type_script
                    build_type_script="$(basename "${bash_script}")"
                    build_type_script="${build_type_script#test}"
                    build_type_script="${build_type_script%.sh}"
                    build_type_script="$(tr '[:upper:]' '[:lower:]' <<< "${build_type_script:0:1}")${build_type_script:1}"
                    CROMWELL_BUILD_TYPE="${build_type_script}"
                    break
                fi
            done
            ;;
    esac

    local backend_type
    backend_type="${CROMWELL_BUILD_TYPE}"
    backend_type="${backend_type#centaurEngineUpgrade}"
    backend_type="${backend_type#centaurWdlUpgrade}"
    backend_type="${backend_type#centaur}"
    backend_type="${backend_type#conformance}"
    backend_type="$(echo "${backend_type}" | sed 's/\([A-Z]\)/_\1/g' | tr '[:upper:]' '[:lower:]' | cut -c 2-)"
    CROMWELL_BUILD_BACKEND_TYPE="${backend_type}"

    CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${CROMWELL_BUILD_BACKEND_TYPE}_application.conf"

    if [ -z "${CROMWELL_BUILD_OPTIONAL_SECURE}" ]; then
        CROMWELL_BUILD_OPTIONAL_SECURE=false
    fi

    if [ -z "${CROMWELL_BUILD_REQUIRES_SECURE}" ]; then
        CROMWELL_BUILD_REQUIRES_SECURE=false
    fi

    if [ -z "${JES_TOKEN}" ]; then
        JES_TOKEN="jes token is not set as an environment variable"
    fi

    CROMWELL_BUILD_RANDOM_256_BITS_BASE64="$(dd bs=1 count=32 if=/dev/urandom 2>/dev/null | base64 | tr -d '\n')"

    local hours_to_minutes
    hours_to_minutes=60
    CROMWELL_BUILD_HEARTBEAT_MINUTES=$((20 * ${hours_to_minutes}))

    export CROMWELL_BUILD_BACKEND_TYPE
    export CROMWELL_BUILD_BRANCH
    export CROMWELL_BUILD_CROMWELL_CONFIG
    export CROMWELL_BUILD_CROMWELL_LOG
    export CROMWELL_BUILD_EVENT
    export CROMWELL_BUILD_EXIT_FUNCTIONS
    export CROMWELL_BUILD_GENERATE_COVERAGE
    export CROMWELL_BUILD_GIT_USER_EMAIL
    export CROMWELL_BUILD_GIT_USER_NAME
    export CROMWELL_BUILD_HEARTBEAT_MESSAGE
    export CROMWELL_BUILD_HEARTBEAT_MINUTES
    export CROMWELL_BUILD_HOME_DIRECTORY
    export CROMWELL_BUILD_IS_CI
    export CROMWELL_BUILD_IS_SECURE
    export CROMWELL_BUILD_LOG_DIRECTORY
    export CROMWELL_BUILD_MYSQL_HOSTNAME
    export CROMWELL_BUILD_MYSQL_PASSWORD
    export CROMWELL_BUILD_MYSQL_PORT
    export CROMWELL_BUILD_MYSQL_SCHEMA
    export CROMWELL_BUILD_MYSQL_USERNAME
    export CROMWELL_BUILD_NUMBER
    export CROMWELL_BUILD_OS
    export CROMWELL_BUILD_OS_DARWIN
    export CROMWELL_BUILD_OS_LINUX
    export CROMWELL_BUILD_PROVIDER
    export CROMWELL_BUILD_PROVIDER_TRAVIS
    export CROMWELL_BUILD_PROVIDER_JENKINS
    export CROMWELL_BUILD_PROVIDER_UNKNOWN
    export CROMWELL_BUILD_RANDOM_256_BITS_BASE64
    export CROMWELL_BUILD_REQUIRES_SECURE
    export CROMWELL_BUILD_OPTIONAL_SECURE
    export CROMWELL_BUILD_RESOURCES_SOURCES
    export CROMWELL_BUILD_RESOURCES_DIRECTORY
    export CROMWELL_BUILD_ROOT_DIRECTORY
    export CROMWELL_BUILD_SCRIPTS_DIRECTORY
    export CROMWELL_BUILD_TAG
    export CROMWELL_BUILD_TYPE
    export CROMWELL_BUILD_URL
    export CROMWELL_BUILD_VAULT_TOKEN
    export CROMWELL_BUILD_WAIT_FOR_IT_BRANCH
    export CROMWELL_BUILD_WAIT_FOR_IT_FILENAME
    export CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT
    export CROMWELL_BUILD_WAIT_FOR_IT_URL
}

cromwell::private::echo_build_variables() {
    echo "CROMWELL_BUILD_IS_CI='${CROMWELL_BUILD_IS_CI}'"
    echo "CROMWELL_BUILD_IS_SECURE='${CROMWELL_BUILD_IS_SECURE}'"
    echo "CROMWELL_BUILD_REQUIRES_SECURE='${CROMWELL_BUILD_REQUIRES_SECURE}'"
    echo "CROMWELL_BUILD_OPTIONAL_SECURE='${CROMWELL_BUILD_OPTIONAL_SECURE}'"
    echo "CROMWELL_BUILD_TYPE='${CROMWELL_BUILD_TYPE}'"
    echo "CROMWELL_BUILD_BRANCH='${CROMWELL_BUILD_BRANCH}'"
    echo "CROMWELL_BUILD_EVENT='${CROMWELL_BUILD_EVENT}'"
    echo "CROMWELL_BUILD_TAG='${CROMWELL_BUILD_TAG}'"
    echo "CROMWELL_BUILD_NUMBER='${CROMWELL_BUILD_NUMBER}'"
    echo "CROMWELL_BUILD_PROVIDER='${CROMWELL_BUILD_PROVIDER}'"
    echo "CROMWELL_BUILD_OS='${CROMWELL_BUILD_OS}'"
    echo "CROMWELL_BUILD_URL='${CROMWELL_BUILD_URL}'"
}

cromwell::private::create_centaur_variables() {
    CROMWELL_BUILD_CENTAUR_TYPE_STANDARD="standard"
    CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION="integration"
    CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE="engineUpgrade"

    if [ -z "${CROMWELL_BUILD_CENTAUR_TYPE}" ]; then
        if [[ "${CROMWELL_BUILD_TYPE}" = centaurEngineUpgrade* ]]; then
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE}"
        else
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_STANDARD}"
        fi
    fi

    CROMWELL_BUILD_CENTAUR_RESOURCES="${CROMWELL_BUILD_ROOT_DIRECTORY}/centaur/src/main/resources"
    CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY="${CROMWELL_BUILD_CENTAUR_RESOURCES}/${CROMWELL_BUILD_CENTAUR_TYPE}TestCases"
    CROMWELL_BUILD_CENTAUR_TEST_RENDERED="${CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY}/rendered"
    CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application.conf"
    CROMWELL_BUILD_CENTAUR_LOG="${CROMWELL_BUILD_LOG_DIRECTORY}/centaur.log"

    case "${CROMWELL_BUILD_CENTAUR_TYPE}" in
        "${CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION}")
            CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=512000
            ;;
        *)
            CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=128000
            ;;
    esac

    export CROMWELL_BUILD_CENTAUR_CONFIG
    export CROMWELL_BUILD_CENTAUR_LOG
    export CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY
    export CROMWELL_BUILD_CENTAUR_TEST_RENDERED
    export CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT
    export CROMWELL_BUILD_CENTAUR_RESOURCES
    export CROMWELL_BUILD_CENTAUR_TYPE
    export CROMWELL_BUILD_CENTAUR_TYPE_STANDARD
    export CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION
    export CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE
}

cromwell::private::verify_secure_build() {
    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
            if [ "${CROMWELL_BUILD_IS_SECURE}" != "true" ] && [ "${CROMWELL_BUILD_REQUIRES_SECURE}" = "true" ]; then
                echo "********************************************************"
                echo "********************************************************"
                echo "**                                                    **"
                echo "**  WARNING: Encrypted keys are unavailable. Exiting. **"
                echo "**                                                    **"
                echo "********************************************************"
                echo "********************************************************"
                exit 0
            fi
            ;;
        *)
            ;;
    esac
}

cromwell::private::export_conformance_variables() {
    CROMWELL_BUILD_CWL_TOOL_VERSION="1.0.20180809224403"
    CROMWELL_BUILD_CWL_TEST_VERSION="1.0.20180601100346"
    CROMWELL_BUILD_CWL_TEST_COMMIT="eb73b5e70e65ab9303a814bd1c230b927018da8f" # use known git hash to avoid changes
    CROMWELL_BUILD_CWL_TEST_RUNNER="${CROMWELL_BUILD_ROOT_DIRECTORY}/centaurCwlRunner/src/bin/centaur-cwl-runner.bash"
    CROMWELL_BUILD_CWL_TEST_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/common-workflow-language"
    CROMWELL_BUILD_CWL_TEST_RESOURCES="${CROMWELL_BUILD_CWL_TEST_DIRECTORY}/v1.0/v1.0"
    CROMWELL_BUILD_CWL_TEST_WDL="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cwl_conformance_test.wdl"
    CROMWELL_BUILD_CWL_TEST_INPUTS="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cwl_conformance_test.inputs.json"
    CROMWELL_BUILD_CWL_TEST_OUTPUT="${CROMWELL_BUILD_LOG_DIRECTORY}/cwl_conformance_test.out.txt"
    CROMWELL_BUILD_CWL_TEST_PARALLELISM=10 # Set too high will cause false negatives due to cromwell server timeouts.

    export CROMWELL_BUILD_CWL_TOOL_VERSION
    export CROMWELL_BUILD_CWL_TEST_VERSION
    export CROMWELL_BUILD_CWL_TEST_COMMIT
    export CROMWELL_BUILD_CWL_TEST_RUNNER
    export CROMWELL_BUILD_CWL_TEST_DIRECTORY
    export CROMWELL_BUILD_CWL_TEST_RESOURCES
    export CROMWELL_BUILD_CWL_TEST_WDL
    export CROMWELL_BUILD_CWL_TEST_INPUTS
    export CROMWELL_BUILD_CWL_TEST_OUTPUT
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

cromwell::private::install_wait_for_it() {
    curl -s "${CROMWELL_BUILD_WAIT_FOR_IT_URL}" > "$CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT"
    chmod +x "$CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT"
}

cromwell::private::create_mysql_cromwell_test() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::create_mysql_cromwell_test
    else
        "$CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT" -t 120 "${CROMWELL_BUILD_MYSQL_HOSTNAME}:${CROMWELL_BUILD_MYSQL_PORT}"
        mysql \
            --host="${CROMWELL_BUILD_MYSQL_HOSTNAME}" \
            --port="${CROMWELL_BUILD_MYSQL_PORT}" \
            --user="${CROMWELL_BUILD_MYSQL_USERNAME}" \
            --password="${CROMWELL_BUILD_MYSQL_PASSWORD}" \
            --execute \
            "SET GLOBAL sql_mode = 'STRICT_ALL_TABLES'; CREATE DATABASE IF NOT EXISTS ${CROMWELL_BUILD_MYSQL_SCHEMA};"
    fi
}

cromwell::private::pull_common_docker_images() {
    # All tests use ubuntu:latest - make sure it's there before starting the tests
    # because pulling the image during some of the tests would cause them to fail
    # (specifically output_redirection which expects a specific value in stderr)
    docker pull ubuntu
}

cromwell::private::install_cwltest() {
    # TODO: No clue why these are needed for cwltool. If you know please update this comment.
    sudo apt-get install procps || true
    sudo -H pip install 'requests[security]'
    sudo -H pip install --ignore-installed cwltool=="${CROMWELL_BUILD_CWL_TOOL_VERSION}"
    sudo -H pip install cwltest=="${CROMWELL_BUILD_CWL_TEST_VERSION}"
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

cromwell::private::write_cwl_test_inputs() {
    cat <<JSON >"${CROMWELL_BUILD_CWL_TEST_INPUTS}"
{
    "cwl_conformance_test.cwl_dir": "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}",
    "cwl_conformance_test.test_result_output": "${CROMWELL_BUILD_CWL_TEST_OUTPUT}",
    "cwl_conformance_test.centaur_cwl_runner": "${CROMWELL_BUILD_CWL_TEST_RUNNER}",
    "cwl_conformance_test.conformance_expected_failures":
        "${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${CROMWELL_BUILD_BACKEND_TYPE}_conformance_expected_failures.txt",
    "cwl_conformance_test.timeout": 1200
}
JSON
}

cromwell::private::docker_login() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::docker_login
    else
        local dockerhub_auth_include
        dockerhub_auth_include="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/dockerhub_auth.inc.sh"
        if [ -f "${dockerhub_auth_include}" ]; then
            # shellcheck source=/dev/null
            source "${dockerhub_auth_include}"
        fi
    fi
}

cromwell::private::vault_login() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::vault_login
    elif [ "${CROMWELL_BUILD_IS_SECURE}" = "true" ]; then
        case "${CROMWELL_BUILD_PROVIDER}" in
            "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
                # Login to vault to access secrets
                local vault_token
                vault_token="${JES_TOKEN}"
                docker run --rm \
                    -v "${CROMWELL_BUILD_HOME_DIRECTORY}:/root:rw" \
                    broadinstitute/dsde-toolbox \
                    vault auth "${vault_token}" < /dev/null > /dev/null && echo vault auth success \
                || true
                unset vault_token
                ;;
            *)
                ;;
        esac
    fi
}

cromwell::private::render_secure_resources() {
    docker run --rm \
        -v "${CROMWELL_BUILD_VAULT_TOKEN}:/root/.vault-token" \
        -v "${CROMWELL_BUILD_RESOURCES_DIRECTORY}:/resources" \
        -e ENVIRONMENT=not_used \
        -e INPUT_PATH=/resources \
        -e OUT_PATH=/resources \
        broadinstitute/dsde-toolbox render-templates.sh \
    || if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
        echo
        echo "Continuing without rendering secure resources."
    else
        echo
        echo "**************************************************************"
        echo "**************************************************************"
        echo "**                                                          **"
        echo "**        WARNING: Unable to render vault resources.        **"
        echo "**  '*.ctmpl' files should be copied and updated manually.  **"
        echo "**                                                          **"
        echo "**************************************************************"
        echo "**************************************************************"
    fi
}

cromwell::private::setup_secure_resources() {
    if [ "${CROMWELL_BUILD_REQUIRES_SECURE}" == "true" ] || [ "${CROMWELL_BUILD_OPTIONAL_SECURE}" == "true" ]; then
        case "${CROMWELL_BUILD_PROVIDER}" in
            "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
                cromwell::private::vault_login
                cromwell::private::render_secure_resources
                cromwell::private::docker_login
                ;;
            "${CROMWELL_BUILD_PROVIDER_JENKINS}")
                ;;
            *)
                cromwell::private::render_secure_resources
                ;;
        esac
    fi
}

cromwell::private::make_build_directories() {
    if [ "${CROMWELL_BUILD_PROVIDER}" == "${CROMWELL_BUILD_PROVIDER_JENKINS}" ]; then
        sudo chmod -R a+w .
    fi
    mkdir -p "${CROMWELL_BUILD_LOG_DIRECTORY}"
    mkdir -p "${CROMWELL_BUILD_RESOURCES_DIRECTORY}"
    cp -r "${CROMWELL_BUILD_RESOURCES_SOURCES}"/* "${CROMWELL_BUILD_RESOURCES_DIRECTORY}"
}

cromwell::private::find_cromwell_jar() {
    CROMWELL_BUILD_CROMWELL_JAR="$( \
        find "${CROMWELL_BUILD_ROOT_DIRECTORY}/server/target/scala-2.12" -name "cromwell-*.jar" \
        | head -n 1 \
        2> /dev/null \
        )"
    export CROMWELL_BUILD_CROMWELL_JAR
}

cromwell::private::setup_prior_version_resources() {
    local current_version
    local prior_version
    local prior_config
    current_version="$( \
        grep 'val cromwellVersion' "${CROMWELL_BUILD_ROOT_DIRECTORY}/project/Version.scala" \
        | awk -F \" '{print $2}' \
        )"
    prior_version=$((current_version - 1))

    CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_JAR="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_${prior_version}.jar"
    export CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_JAR

    prior_config="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${CROMWELL_BUILD_BACKEND_TYPE}_${prior_version}_application.conf"
    if [ -f "${prior_config}" ]; then
        CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_CONFIG="${prior_config}"
        export CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_CONFIG
    fi

    docker run \
      --rm \
      --entrypoint= \
      -v "${PWD}:${PWD}" \
      broadinstitute/cromwell:"${prior_version}" \
      cp /app/cromwell.jar "${CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_JAR}"
}

cromwell::private::exists_cromwell_jar() {
    test -s "${CROMWELL_BUILD_CROMWELL_JAR}"
}

cromwell::private::assemble_jars() {
    CROMWELL_SBT_ASSEMBLY_LOG_LEVEL=error sbt coverage assembly -error
}

cromwell::private::generate_code_coverage() {
    sbt coverageReport -warn
    sbt coverageAggregate -warn
    bash <(curl -s https://codecov.io/bash) > /dev/null || true
}

cromwell::private::publish_artifacts_only() {
    CROMWELL_SBT_ASSEMBLY_LOG_LEVEL=warn sbt "$@" publish -warn
}

cromwell::private::publish_artifacts_and_docker() {
    CROMWELL_SBT_ASSEMBLY_LOG_LEVEL=warn sbt "$@" publish dockerBuildAndPush -warn
}

cromwell::private::publish_artifacts_check() {
    sbt verifyArtifactoryCredentialsExist -warn
}

# Some CI environments want to know when new docker images are published. They do not currently poll dockerhub but do
# poll github. To help those environments, signal that a new set of docker images has been published to dockerhub by
# updating a well known branch in github.
cromwell::private::push_publish_complete() {
    local github_private_deploy_key="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/github_private_deploy_key"
    local git_repo="git@github.com:broadinstitute/cromwell.git"
    local git_publish_branch="${CROMWELL_BUILD_BRANCH}_publish_complete"
    local git_publish_remote="publish_complete"
    local git_publish_message="publish complete [skip ci]"

    # Loosely adapted from https://github.com/broadinstitute/workbench-libs/blob/435a932/scripts/version_update.sh
    mkdir publish_complete
    (
        cd publish_complete || exit 1

        git init
        git config core.sshCommand "ssh -i ${github_private_deploy_key} -F /dev/null"
        git config user.email "${CROMWELL_BUILD_GIT_USER_EMAIL}"
        git config user.name "${CROMWELL_BUILD_GIT_USER_NAME}"

        git remote add "${git_publish_remote}" "${git_repo}"
        git checkout -b "${git_publish_branch}"
        git commit --allow-empty -m "${git_publish_message}"
        git push -f "${git_publish_remote}" "${git_publish_branch}"
    )
}

cromwell::private::start_build_heartbeat() {
    # Sleep one minute between printouts, but don't zombie forever
    for ((i=0; i < "${CROMWELL_BUILD_HEARTBEAT_MINUTES}"; i++)); do
        sleep 60
        printf "${CROMWELL_BUILD_HEARTBEAT_MESSAGE}"
    done &
    CROMWELL_BUILD_HEARTBEAT_PID=$!
}

cromwell::private::start_cromwell_log_tail() {
    while [ ! -f "${CROMWELL_BUILD_CROMWELL_LOG}" ]; do
        sleep 2
    done && tail -n 0 -f "${CROMWELL_BUILD_CROMWELL_LOG}" 2> /dev/null &
    CROMWELL_BUILD_CROMWELL_LOG_TAIL_PID=$!
}

cromwell::private::start_centaur_log_tail() {
    while [ ! -f "${CROMWELL_BUILD_CENTAUR_LOG}" ]; do
        sleep 2
    done && tail -n 0 -f "${CROMWELL_BUILD_CENTAUR_LOG}" 2> /dev/null &
    CROMWELL_BUILD_CENTAUR_LOG_TAIL_PID=$!
}

cromwell::private::cat_centaur_log() {
    echo "CENTAUR LOG"
    cat "${CROMWELL_BUILD_CENTAUR_LOG}"
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
    cromwell::private::exec_test_script
}

cromwell::build::setup_common_environment() {
    cromwell::private::check_debug
    cromwell::private::create_build_variables
    cromwell::private::echo_build_variables
    cromwell::private::verify_secure_build
    cromwell::private::make_build_directories
    cromwell::private::setup_secure_resources

    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
            cromwell::private::delete_boto_config
            cromwell::private::delete_sbt_boot
            cromwell::private::upgrade_pip
            cromwell::private::pull_common_docker_images
            cromwell::private::install_wait_for_it
            cromwell::private::create_mysql_cromwell_test
            ;;
        "${CROMWELL_BUILD_PROVIDER_JENKINS}")
            cromwell::private::delete_boto_config
            cromwell::private::delete_sbt_boot
            cromwell::private::upgrade_pip
            cromwell::private::install_wait_for_it
            cromwell::private::create_mysql_cromwell_test
            ;;
        *)
            cromwell::private::pull_common_docker_images
            ;;
    esac
}

cromwell::build::setup_centaur_environment() {
    cromwell::private::create_centaur_variables
    if [ "${CROMWELL_BUILD_CENTAUR_TYPE}" = "${CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE}" ]; then
        cromwell::private::setup_prior_version_resources
    fi
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

cromwell::build::run_centaur() {
    "${CROMWELL_BUILD_ROOT_DIRECTORY}/centaur/test_cromwell.sh" \
        -n "${CROMWELL_BUILD_CENTAUR_CONFIG}" \
        -l "${CROMWELL_BUILD_LOG_DIRECTORY}" \
        -g \
        "$@"
}

cromwell::build::setup_conformance_environment() {
    cromwell::private::export_conformance_variables
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ]; then
        cromwell::private::install_cwltest
    fi
    cromwell::private::checkout_pinned_cwl
    cromwell::private::write_cwl_test_inputs
    cromwell::private::start_build_heartbeat
    cromwell::private::add_exit_function cromwell::private::cat_conformance_log
    cromwell::private::add_exit_function cromwell::private::kill_build_heartbeat
}

cromwell::build::assemble_jars() {
    cromwell::private::find_cromwell_jar
    if [ "${CROMWELL_BUILD_IS_CI}" = "true" ] || ! cromwell::private::exists_cromwell_jar; then
        echo "Please wait, building jars…"
        cromwell::private::assemble_jars
    fi
    cromwell::private::find_cromwell_jar
    if ! cromwell::private::exists_cromwell_jar; then
        echo "Error: find_cromwell_jar did not locate a cromwell jar even after assembly" >&2
        exit 1
    fi
}

cromwell::build::generate_code_coverage() {
    if [ "${CROMWELL_BUILD_GENERATE_COVERAGE}" = "true" ]; then
        cromwell::private::generate_code_coverage
    fi
}

cromwell::build::publish_artifacts() {
    if [ "${CROMWELL_BUILD_PROVIDER}" == "${CROMWELL_BUILD_PROVIDER_TRAVIS}" ] && \
        [ "${CROMWELL_BUILD_TYPE}" == "sbt" ] && \
        [ "${CROMWELL_BUILD_EVENT}" == "push" ]; then

        if [ "${CROMWELL_BUILD_BRANCH}" = "develop" ]; then
            # Publish images for both the "cromwell develop branch" and the "cromwell dev environment".
            CROMWELL_SBT_DOCKER_TAGS=develop,dev \
                cromwell::private::publish_artifacts_and_docker \
                -Dproject.isSnapshot=true
            cromwell::private::push_publish_complete

        elif [[ "${CROMWELL_BUILD_BRANCH}" =~ ^[0-9\.]+_hotfix$ ]]; then
            # Docker tags float. "30" is the latest hotfix. Those dockers are published here on each hotfix commit.
            cromwell::private::publish_artifacts_and_docker -Dproject.isSnapshot=false

        elif [ -n "${CROMWELL_BUILD_TAG:+set}" ]; then
            # Artifact tags are static. Once "30" is set that is only "30" forever. Those artifacts are published here.
            cromwell::private::publish_artifacts_only \
                -Dproject.version="${CROMWELL_BUILD_TAG}" \
                -Dproject.isSnapshot=false

        elif [ "${CROMWELL_BUILD_IS_SECURE}" = "true" ]; then
            cromwell::private::publish_artifacts_check

        fi

    fi
}

cromwell::build::add_exit_function() {
    cromwell::private::add_exit_function "$1"
}

cromwell::build::kill_tree() {
    cromwell::private::kill_tree "$1"
}
