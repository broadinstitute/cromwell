#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

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
#     Quick debug scripts. Example: `crmdbg=y src/ci/bin/testCentaurLocal.sh`
#
#   - crmcit
#     Simulate a centaur integration test build. Example: `crmcit=y src/ci/bin/testCentaurPapiV2beta.sh`
#
#   - crmddm
#     Use "Docker Desktop for Mac" DNS names instead of `localhost`.
#     Example: `crmddm=y src/ci/bin/testCentaurHoricromtalPapiV2alpha1.sh`
#     More info: https://docs.docker.com/docker-for-mac/networking/#i-want-to-connect-from-a-container-to-a-service-on-the-host

cromwell::private::check_debug() {
    # shellcheck disable=SC2154
    if [[ -n "${crmdbg:+set}" ]]; then
        set -o xtrace
    fi

    # shellcheck disable=SC2154
    if [[ -n "${crmcit:+set}" ]]; then
        CROMWELL_BUILD_CENTAUR_TYPE="integration"
    fi

    # shellcheck disable=SC2154
    if [[ -n "${crmddm:+set}" ]]; then
        CROMWELL_BUILD_DOCKER_LOCALHOST="host.docker.internal"
    fi
}

# Exports environment variables used for scripts.
cromwell::private::create_build_variables() {
    CROMWELL_BUILD_PROVIDER_GITHUB="github"
    CROMWELL_BUILD_PROVIDER_UNKNOWN="unknown"

    # CI providers set an env like `GITHUB_ACTIONS` that indicates to jobs where they're running
    # https://docs.github.com/en/actions/learn-github-actions/variables#default-environment-variables
    if [[ "${GITHUB_ACTIONS-false}" == "true" ]]; then
        CROMWELL_BUILD_PROVIDER="${CROMWELL_BUILD_PROVIDER_GITHUB}"
    else
        CROMWELL_BUILD_PROVIDER="${CROMWELL_BUILD_PROVIDER_UNKNOWN}"
    fi

    # simplified from https://stackoverflow.com/a/18434831/3320205
    CROMWELL_BUILD_OS_DARWIN="darwin";
    CROMWELL_BUILD_OS_LINUX="linux";
    case "${OSTYPE-unknown}" in
        darwin*)  CROMWELL_BUILD_OS="${CROMWELL_BUILD_OS_DARWIN}" ;;
        linux*)   CROMWELL_BUILD_OS="${CROMWELL_BUILD_OS_LINUX}" ;;
        *)        CROMWELL_BUILD_OS="unknown_os" ;;
    esac

    CROMWELL_BUILD_HOME_DIRECTORY="${HOME}"
    CROMWELL_BUILD_ROOT_DIRECTORY="$(pwd)"
    CROMWELL_BUILD_LOG_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/target/ci/logs"
    CROMWELL_BUILD_CROMWELL_LOG="${CROMWELL_BUILD_LOG_DIRECTORY}/cromwell.log"

    CROMWELL_BUILD_DOCKER_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/src/ci/docker-compose"
    CROMWELL_BUILD_SCRIPTS_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/src/ci/bin"
    CROMWELL_BUILD_RESOURCES_SOURCES="${CROMWELL_BUILD_ROOT_DIRECTORY}/src/ci/resources"
    CROMWELL_BUILD_RESOURCES_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/target/ci/resources"

    CROMWELL_BUILD_GIT_SECRETS_DIRECTORY="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/git-secrets"
    CROMWELL_BUILD_GIT_SECRETS_COMMIT="ad82d68ee924906a0401dfd48de5057731a9bc84"
    CROMWELL_BUILD_WAIT_FOR_IT_FILENAME="wait-for-it.sh"
    CROMWELL_BUILD_WAIT_FOR_IT_BRANCH="db049716e42767d39961e95dd9696103dca813f1"
    CROMWELL_BUILD_WAIT_FOR_IT_URL="https://raw.githubusercontent.com/vishnubob/wait-for-it/${CROMWELL_BUILD_WAIT_FOR_IT_BRANCH}/${CROMWELL_BUILD_WAIT_FOR_IT_FILENAME}"
    CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${CROMWELL_BUILD_WAIT_FOR_IT_FILENAME}"
    CROMWELL_BUILD_VAULT_ZIP="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/vault.zip"
    CROMWELL_BUILD_VAULT_EXECUTABLE="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/vault"
    CROMWELL_BUILD_EXIT_FUNCTIONS="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_build_exit_functions.$$"

    if [[ -n "${VIRTUAL_ENV:+set}" ]]; then
        CROMWELL_BUILD_IS_VIRTUAL_ENV=true
    else
        CROMWELL_BUILD_IS_VIRTUAL_ENV=false
    fi

    CROMWELL_BUILD_CURRENT_VERSION_NUMBER="$( \
        grep 'val cromwellVersion' "${CROMWELL_BUILD_ROOT_DIRECTORY}/project/Version.scala" \
        | awk -F \" '{print $2}' \
        )"

    if git merge-base --is-ancestor "${CROMWELL_BUILD_CURRENT_VERSION_NUMBER}" HEAD 2>/dev/null; then
        CROMWELL_BUILD_IS_HOTFIX=true
    else
        CROMWELL_BUILD_IS_HOTFIX=false
    fi

    if [[ "${CROMWELL_BUILD_IS_HOTFIX}" == "true" ]]; then
        CROMWELL_BUILD_PRIOR_VERSION_NUMBER=${CROMWELL_BUILD_CURRENT_VERSION_NUMBER}
    else
        CROMWELL_BUILD_PRIOR_VERSION_NUMBER=$((CROMWELL_BUILD_CURRENT_VERSION_NUMBER - 1))
    fi

    local git_revision
    if git_revision="$(git rev-parse --short=7 HEAD 2>/dev/null)"; then
        CROMWELL_BUILD_GIT_HASH_SUFFIX="g${git_revision}"
    else
        CROMWELL_BUILD_GIT_HASH_SUFFIX="gUNKNOWN"
    fi

    # PR #6790 disabled the conditional that skips tests for documentation-only PRs, because
    # those PRs (and only those PRs) were uniformly failing tests with a nondescript error.
    # https://broadinstitute.slack.com/archives/GHYJZ2ZE0/p1656625952888149?thread_ts=1656620572.975059&cid=GHYJZ2ZE0

    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_GITHUB}")
            CROMWELL_BUILD_IS_CI=true
            CROMWELL_BUILD_TYPE="${BUILD_TYPE}"
            CROMWELL_BUILD_BRANCH="${GITHUB_REF_NAME}"
            CROMWELL_BUILD_EVENT="${GITHUB_EVENT_NAME}"
            CROMWELL_BUILD_TAG=""
            CROMWELL_BUILD_NUMBER="${GITHUB_RUN_ID}"
            CROMWELL_BUILD_URL="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${GITHUB_RUN_ID}"
            CROMWELL_BUILD_GIT_USER_EMAIL=""
            CROMWELL_BUILD_GIT_USER_NAME="${GITHUB_ACTOR}"
            CROMWELL_BUILD_HEARTBEAT_PATTERN="…"
            CROMWELL_BUILD_GENERATE_COVERAGE=true
            CROMWELL_BUILD_RUN_TESTS=true
            ;;
        *)
            CROMWELL_BUILD_IS_CI=false
            CROMWELL_BUILD_TYPE="unknown"
            CROMWELL_BUILD_BRANCH="unknown"
            CROMWELL_BUILD_EVENT="unknown"
            CROMWELL_BUILD_TAG=""
            CROMWELL_BUILD_NUMBER=""
            CROMWELL_BUILD_URL=""
            CROMWELL_BUILD_GIT_USER_EMAIL="unknown.git.user@example.org"
            CROMWELL_BUILD_GIT_USER_NAME="Unknown Git User"
            CROMWELL_BUILD_HEARTBEAT_PATTERN="…"
            CROMWELL_BUILD_GENERATE_COVERAGE="${CROMWELL_BUILD_GENERATE_COVERAGE:-true}"
            CROMWELL_BUILD_RUN_TESTS=true

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
    backend_type="${backend_type#centaurPapiUpgrade}"
    backend_type="${backend_type#centaurWdlUpgrade}"
    backend_type="${backend_type#centaurHoricromtal}"
    backend_type="${backend_type#centaur}"
    backend_type="$(echo "${backend_type}" | sed 's/\([A-Z]\)/_\1/g' | tr '[:upper:]' '[:lower:]' | cut -c 2-)"
    CROMWELL_BUILD_BACKEND_TYPE="${backend_type}"

    # It may be possible to trim this down to e.g. `server/assembly` for some jobs
    CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND="assembly"

    if [[ "${CROMWELL_BUILD_GENERATE_COVERAGE}" == "true" ]]; then
        CROMWELL_BUILD_SBT_COVERAGE_COMMAND="coverage"
    else
        CROMWELL_BUILD_SBT_COVERAGE_COMMAND=""
    fi

    CROMWELL_BUILD_SBT_INCLUDE="${BUILD_SBT_INCLUDE:-}"
    CROMWELL_BUILD_SBT_EXCLUDE="${BUILD_SBT_EXCLUDE:-}"

    case "${CROMWELL_BUILD_TYPE}" in
        centaurHoricromtalPapiV2beta*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2beta_horicromtal_application.conf"
            ;;
        centaurBlob*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_blob_test.conf"
            ;;
        *)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${CROMWELL_BUILD_BACKEND_TYPE}_application.conf"
            ;;
    esac

    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        CROMWELL_BUILD_DOCKER_TAG="${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_NUMBER}"
    else
        CROMWELL_BUILD_DOCKER_TAG="${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_TYPE}-${CROMWELL_BUILD_GIT_HASH_SUFFIX}"
    fi

    # Trim and replace invalid characters in the docker tag
    # https://docs.docker.com/engine/reference/commandline/tag/#extended-description
    CROMWELL_BUILD_DOCKER_TAG="${CROMWELL_BUILD_DOCKER_TAG:0:128}"
    CROMWELL_BUILD_DOCKER_TAG="${CROMWELL_BUILD_DOCKER_TAG//[^a-zA-Z0-9.-]/_}"

    CROMWELL_BUILD_REQUIRES_PRIOR_VERSION="${CROMWELL_BUILD_REQUIRES_PRIOR_VERSION-false}"

    local hours_to_minutes
    hours_to_minutes=60
    CROMWELL_BUILD_HEARTBEAT_MINUTES=$((20 * hours_to_minutes))

    export CROMWELL_BUILD_BACKEND_TYPE
    export CROMWELL_BUILD_BRANCH
    export CROMWELL_BUILD_CROMWELL_CONFIG
    export CROMWELL_BUILD_CROMWELL_LOG
    export CROMWELL_BUILD_CURRENT_VERSION_NUMBER
    export CROMWELL_BUILD_DOCKER_DIRECTORY
    export CROMWELL_BUILD_DOCKER_TAG
    export CROMWELL_BUILD_EVENT
    export CROMWELL_BUILD_EXIT_FUNCTIONS
    export CROMWELL_BUILD_GENERATE_COVERAGE
    export CROMWELL_BUILD_GIT_HASH_SUFFIX
    export CROMWELL_BUILD_GIT_SECRETS_COMMIT
    export CROMWELL_BUILD_GIT_SECRETS_DIRECTORY
    export CROMWELL_BUILD_GIT_USER_EMAIL
    export CROMWELL_BUILD_GIT_USER_NAME
    export CROMWELL_BUILD_HEARTBEAT_MINUTES
    export CROMWELL_BUILD_HEARTBEAT_PATTERN
    export CROMWELL_BUILD_HOME_DIRECTORY
    export CROMWELL_BUILD_IS_CI
    export CROMWELL_BUILD_IS_HOTFIX
    export CROMWELL_BUILD_IS_VIRTUAL_ENV
    export CROMWELL_BUILD_LOG_DIRECTORY
    export CROMWELL_BUILD_NUMBER
    export CROMWELL_BUILD_OS
    export CROMWELL_BUILD_OS_DARWIN
    export CROMWELL_BUILD_OS_LINUX
    export CROMWELL_BUILD_PRIOR_VERSION_NUMBER
    export CROMWELL_BUILD_PROVIDER
    export CROMWELL_BUILD_PROVIDER_UNKNOWN
    export CROMWELL_BUILD_REQUIRES_PRIOR_VERSION
    export CROMWELL_BUILD_REQUIRES_SECURE
    export CROMWELL_BUILD_RESOURCES_DIRECTORY
    export CROMWELL_BUILD_RESOURCES_SOURCES
    export CROMWELL_BUILD_ROOT_DIRECTORY
    export CROMWELL_BUILD_RUN_TESTS
    export CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND
    export CROMWELL_BUILD_SBT_COVERAGE_COMMAND
    export CROMWELL_BUILD_SBT_EXCLUDE
    export CROMWELL_BUILD_SBT_INCLUDE
    export CROMWELL_BUILD_SCRIPTS_DIRECTORY
    export CROMWELL_BUILD_TAG
    export CROMWELL_BUILD_TYPE
    export CROMWELL_BUILD_UNIT_TEST_EXCLUDE_TAGS
    export CROMWELL_BUILD_URL
    export CROMWELL_BUILD_VAULT_EXECUTABLE
    export CROMWELL_BUILD_VAULT_ZIP
    export CROMWELL_BUILD_WAIT_FOR_IT_BRANCH
    export CROMWELL_BUILD_WAIT_FOR_IT_FILENAME
    export CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT
    export CROMWELL_BUILD_WAIT_FOR_IT_URL
}

cromwell::private::echo_build_variables() {
    echo "CROMWELL_BUILD_IS_CI='${CROMWELL_BUILD_IS_CI}'"
    echo "CROMWELL_BUILD_TYPE='${CROMWELL_BUILD_TYPE}'"
    echo "CROMWELL_BUILD_BRANCH='${CROMWELL_BUILD_BRANCH}'"
    echo "CROMWELL_BUILD_IS_HOTFIX='${CROMWELL_BUILD_IS_HOTFIX}'"
    echo "CROMWELL_BUILD_CURRENT_VERSION_NUMBER='${CROMWELL_BUILD_CURRENT_VERSION_NUMBER}'"
    echo "CROMWELL_BUILD_PRIOR_VERSION_NUMBER='${CROMWELL_BUILD_PRIOR_VERSION_NUMBER}'"
    echo "CROMWELL_BUILD_EVENT='${CROMWELL_BUILD_EVENT}'"
    echo "CROMWELL_BUILD_TAG='${CROMWELL_BUILD_TAG}'"
    echo "CROMWELL_BUILD_NUMBER='${CROMWELL_BUILD_NUMBER}'"
    echo "CROMWELL_BUILD_PROVIDER='${CROMWELL_BUILD_PROVIDER}'"
    echo "CROMWELL_BUILD_OS='${CROMWELL_BUILD_OS}'"
    echo "CROMWELL_BUILD_URL='${CROMWELL_BUILD_URL}'"
}

# Create environment variables used by the DatabaseTestKit and cromwell::private::create_centaur_variables()
cromwell::private::create_database_variables() {
    CROMWELL_BUILD_DATABASE_USERNAME="cromwell"
    CROMWELL_BUILD_DATABASE_PASSWORD="test"
    CROMWELL_BUILD_DATABASE_SCHEMA="cromwell_test"

    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_GITHUB}")
            CROMWELL_BUILD_MARIADB_HOSTNAME="localhost"
            CROMWELL_BUILD_MARIADB_PORT="23306"
            CROMWELL_BUILD_MARIADB_DOCKER_TAG="${BUILD_MARIADB-}"
            CROMWELL_BUILD_MARIADB_LATEST_HOSTNAME="localhost"
            CROMWELL_BUILD_MARIADB_LATEST_PORT="33306"
            CROMWELL_BUILD_MARIADB_LATEST_TAG="${BUILD_MARIADB_LATEST-}"
            CROMWELL_BUILD_MYSQL_HOSTNAME="localhost"
            CROMWELL_BUILD_MYSQL_PORT="3306"
            CROMWELL_BUILD_MYSQL_DOCKER_TAG="${BUILD_MYSQL-}"
            CROMWELL_BUILD_MYSQL_LATEST_HOSTNAME="localhost"
            CROMWELL_BUILD_MYSQL_LATEST_PORT="13306"
            CROMWELL_BUILD_MYSQL_LATEST_TAG="${BUILD_MYSQL_LATEST-}"
            CROMWELL_BUILD_POSTGRESQL_HOSTNAME="localhost"
            CROMWELL_BUILD_POSTGRESQL_PORT="5432"
            CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG="${BUILD_POSTGRESQL-}"
            CROMWELL_BUILD_POSTGRESQL_LATEST_HOSTNAME="localhost"
            CROMWELL_BUILD_POSTGRESQL_LATEST_PORT="15432"
            CROMWELL_BUILD_POSTGRESQL_LATEST_TAG="${BUILD_POSTGRESQL_LATEST-}"
            ;;
        *)
            if [[ -z "${CROMWELL_BUILD_DOCKER_LOCALHOST-}" ]]; then
                CROMWELL_BUILD_DOCKER_LOCALHOST="localhost"
            fi

            CROMWELL_BUILD_MARIADB_HOSTNAME="${CROMWELL_BUILD_MARIADB_HOSTNAME-${CROMWELL_BUILD_DOCKER_LOCALHOST}}"
            CROMWELL_BUILD_MARIADB_PORT="${CROMWELL_BUILD_MARIADB_PORT-13306}"
            CROMWELL_BUILD_MARIADB_DOCKER_TAG=""
            CROMWELL_BUILD_MARIADB_LATEST_HOSTNAME="${CROMWELL_BUILD_MARIADB_LATEST_HOSTNAME-${CROMWELL_BUILD_DOCKER_LOCALHOST}}"
            CROMWELL_BUILD_MARIADB_LATEST_PORT="${CROMWELL_BUILD_MARIADB_LATEST_PORT-13306}"
            CROMWELL_BUILD_MARIADB_LATEST_TAG=""
            CROMWELL_BUILD_MYSQL_HOSTNAME="${CROMWELL_BUILD_MYSQL_HOSTNAME-${CROMWELL_BUILD_DOCKER_LOCALHOST}}"
            CROMWELL_BUILD_MYSQL_PORT="${CROMWELL_BUILD_MYSQL_PORT-3306}"
            CROMWELL_BUILD_MYSQL_DOCKER_TAG=""
            CROMWELL_BUILD_MYSQL_LATEST_HOSTNAME="${CROMWELL_BUILD_MYSQL_LATEST_HOSTNAME-${CROMWELL_BUILD_DOCKER_LOCALHOST}}"
            CROMWELL_BUILD_MYSQL_LATEST_PORT="${CROMWELL_BUILD_MYSQL_LATEST_PORT-13306}"
            CROMWELL_BUILD_MYSQL_LATEST_TAG=""
            CROMWELL_BUILD_POSTGRESQL_HOSTNAME="${CROMWELL_BUILD_POSTGRESQL_HOSTNAME-${CROMWELL_BUILD_DOCKER_LOCALHOST}}"
            CROMWELL_BUILD_POSTGRESQL_PORT="${CROMWELL_BUILD_POSTGRESQL_PORT-5432}"
            CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG=""
            CROMWELL_BUILD_POSTGRESQL_LATEST_HOSTNAME="${CROMWELL_BUILD_POSTGRESQL_LATEST_HOSTNAME-${CROMWELL_BUILD_DOCKER_LOCALHOST}}"
            CROMWELL_BUILD_POSTGRESQL_LATEST_PORT="${CROMWELL_BUILD_POSTGRESQL_LATEST_PORT-13306}"
            CROMWELL_BUILD_POSTGRESQL_LATEST_TAG=""
            ;;
    esac

    export CROMWELL_BUILD_DATABASE_USERNAME
    export CROMWELL_BUILD_DATABASE_PASSWORD
    export CROMWELL_BUILD_DATABASE_SCHEMA
    export CROMWELL_BUILD_MARIADB_DOCKER_TAG
    export CROMWELL_BUILD_MARIADB_HOSTNAME
    export CROMWELL_BUILD_MARIADB_LATEST_HOSTNAME
    export CROMWELL_BUILD_MARIADB_LATEST_PORT
    export CROMWELL_BUILD_MARIADB_LATEST_TAG
    export CROMWELL_BUILD_MARIADB_PORT
    export CROMWELL_BUILD_MYSQL_DOCKER_TAG
    export CROMWELL_BUILD_MYSQL_HOSTNAME
    export CROMWELL_BUILD_MYSQL_LATEST_HOSTNAME
    export CROMWELL_BUILD_MYSQL_LATEST_PORT
    export CROMWELL_BUILD_MYSQL_LATEST_TAG
    export CROMWELL_BUILD_MYSQL_PORT
    export CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG
    export CROMWELL_BUILD_POSTGRESQL_HOSTNAME
    export CROMWELL_BUILD_POSTGRESQL_LATEST_HOSTNAME
    export CROMWELL_BUILD_POSTGRESQL_LATEST_PORT
    export CROMWELL_BUILD_POSTGRESQL_LATEST_TAG
    export CROMWELL_BUILD_POSTGRESQL_PORT
}

cromwell::private::create_centaur_variables() {
    CROMWELL_BUILD_CENTAUR_TYPE_STANDARD="standard"
    CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION="integration"
    CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE="engineUpgrade"
    CROMWELL_BUILD_CENTAUR_TYPE_PAPI_UPGRADE="papiUpgrade"
    CROMWELL_BUILD_CENTAUR_TYPE_PAPI_UPGRADE_NEW_WORKFLOWS="papiUpgradeNewWorkflows"
    CROMWELL_BUILD_CENTAUR_TYPE_HORICROMTAL_ENGINE_UPGRADE="horicromtalEngineUpgrade"
    CROMWELL_BUILD_CENTAUR_TYPE_HORICROMTAL="horicromtal"
    CROMWELL_BUILD_CENTAUR_TYPE_AZURE_BLOB="azureBlob"

    case "${CROMWELL_BUILD_TYPE}" in
        centaurEngineUpgrade*)
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE}"
            ;;
        centaurPapiUpgradeNewWorkflows*)
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_PAPI_UPGRADE_NEW_WORKFLOWS}"
            ;;
        centaurPapiUpgrade*)
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_PAPI_UPGRADE}"
            ;;
        centaurHoricromtalEngineUpgrade*)
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_HORICROMTAL_ENGINE_UPGRADE}"
            ;;
        centaurHoricromtal*)
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_HORICROMTAL}"
            ;;
        centaurBlob*)
            CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_AZURE_BLOB}"
            ;;
        *)
            # Only set the type if Jenkins, etc. has not already set the centaur type
            if [[ -z "${CROMWELL_BUILD_CENTAUR_TYPE-}" ]]; then
                CROMWELL_BUILD_CENTAUR_TYPE="${CROMWELL_BUILD_CENTAUR_TYPE_STANDARD}"
            fi
            ;;
    esac

    CROMWELL_BUILD_CENTAUR_RESOURCES="${CROMWELL_BUILD_ROOT_DIRECTORY}/centaur/src/main/resources"
    case "${CROMWELL_BUILD_CENTAUR_TYPE}" in
        "${CROMWELL_BUILD_CENTAUR_TYPE_HORICROMTAL}")
            # Use the standard test cases despite the horicromtal Centaur build type.
            CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY="${CROMWELL_BUILD_CENTAUR_RESOURCES}/standardTestCases"
            CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application_horicromtal.conf"
            ;;
        "${CROMWELL_BUILD_CENTAUR_TYPE_HORICROMTAL_ENGINE_UPGRADE}")
            # Use the engine upgrade test cases despite the horicromtal Centaur build type.
            CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY="${CROMWELL_BUILD_CENTAUR_RESOURCES}/engineUpgradeTestCases"
            # Special horicromtal engine upgrade Centaur config with horicromtal assertions turned off.
            CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application_horicromtal_no_assert.conf"
            ;;
        *)
            CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY="${CROMWELL_BUILD_CENTAUR_RESOURCES}/${CROMWELL_BUILD_CENTAUR_TYPE}TestCases"
            CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application.conf"
            ;;
    esac

    CROMWELL_BUILD_CENTAUR_LOG="${CROMWELL_BUILD_LOG_DIRECTORY}/centaur.log"

    local mariadb_jdbc_url
    local mysql_jdbc_url
    local postgresql_jdbc_url

    mariadb_jdbc_url="jdbc:mariadb://${CROMWELL_BUILD_MARIADB_HOSTNAME}:${CROMWELL_BUILD_MARIADB_PORT}/${CROMWELL_BUILD_DATABASE_SCHEMA}?rewriteBatchedStatements=true"
    mysql_jdbc_url="jdbc:mysql://${CROMWELL_BUILD_MYSQL_HOSTNAME}:${CROMWELL_BUILD_MYSQL_PORT}/${CROMWELL_BUILD_DATABASE_SCHEMA}?allowPublicKeyRetrieval=true&useSSL=false&rewriteBatchedStatements=true&serverTimezone=UTC&useInformationSchema=true"
    postgresql_jdbc_url="jdbc:postgresql://${CROMWELL_BUILD_POSTGRESQL_HOSTNAME}:${CROMWELL_BUILD_POSTGRESQL_PORT}/${CROMWELL_BUILD_DATABASE_SCHEMA}?reWriteBatchedInserts=true"

    # Pick **one** of the databases to run Centaur against
    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_GITHUB}")

            if [[ -n "${CROMWELL_BUILD_MYSQL_DOCKER_TAG:+set}" ]]; then
                CROMWELL_BUILD_CENTAUR_SLICK_PROFILE="slick.jdbc.MySQLProfile$"
                CROMWELL_BUILD_CENTAUR_JDBC_DRIVER="com.mysql.cj.jdbc.Driver"
                CROMWELL_BUILD_CENTAUR_JDBC_URL="${mysql_jdbc_url}"

            elif [[ -n "${CROMWELL_BUILD_MARIADB_DOCKER_TAG:+set}" ]]; then
                CROMWELL_BUILD_CENTAUR_SLICK_PROFILE="slick.jdbc.MySQLProfile$"
                CROMWELL_BUILD_CENTAUR_JDBC_DRIVER="org.mariadb.jdbc.Driver"
                CROMWELL_BUILD_CENTAUR_JDBC_URL="${mariadb_jdbc_url}"

            elif [[ -n "${CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG:+set}" ]]; then
                CROMWELL_BUILD_CENTAUR_SLICK_PROFILE="slick.jdbc.PostgresProfile$"
                CROMWELL_BUILD_CENTAUR_JDBC_DRIVER="org.postgresql.Driver"
                CROMWELL_BUILD_CENTAUR_JDBC_URL="${postgresql_jdbc_url}"

            else
                echo "Error: Unable to determine which RDBMS to use for Centaur." >&2
                exit 1

            fi

            CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS=
            ;;
        *)
            CROMWELL_BUILD_CENTAUR_SLICK_PROFILE="${CROMWELL_BUILD_CENTAUR_SLICK_PROFILE-slick.jdbc.MySQLProfile\$}"
            CROMWELL_BUILD_CENTAUR_JDBC_DRIVER="${CROMWELL_BUILD_CENTAUR_JDBC_DRIVER-com.mysql.cj.jdbc.Driver}"
            CROMWELL_BUILD_CENTAUR_JDBC_URL="${CROMWELL_BUILD_CENTAUR_JDBC_URL-${mysql_jdbc_url}}"
            CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS=
            ;;
    esac

    case "${CROMWELL_BUILD_CENTAUR_TYPE}" in
        "${CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION}")
            CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=512000
            CROMWELL_BUILD_CENTAUR_MAX_WORKFLOW_LENGTH="10 hours"
            ;;
        *)
            CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=128000
            CROMWELL_BUILD_CENTAUR_MAX_WORKFLOW_LENGTH="90 minutes"
            ;;
    esac

    CROMWELL_BUILD_CENTAUR_PRIOR_SLICK_PROFILE="${CROMWELL_BUILD_CENTAUR_PRIOR_SLICK_PROFILE-${CROMWELL_BUILD_CENTAUR_SLICK_PROFILE}}"
    CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_DRIVER="${CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_DRIVER-${CROMWELL_BUILD_CENTAUR_JDBC_DRIVER}}"
    CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_URL="${CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_URL-${CROMWELL_BUILD_CENTAUR_JDBC_URL}}"

    CROMWELL_BUILD_CENTAUR_256_BITS_KEY="$(dd bs=1 count=32 if=/dev/urandom 2> /dev/null | base64 | tr -d '\n')"

    export CROMWELL_BUILD_CENTAUR_256_BITS_KEY
    export CROMWELL_BUILD_CENTAUR_CONFIG
    export CROMWELL_BUILD_CENTAUR_JDBC_DRIVER
    export CROMWELL_BUILD_CENTAUR_JDBC_URL
    export CROMWELL_BUILD_CENTAUR_LOG
    export CROMWELL_BUILD_CENTAUR_MAX_WORKFLOW_LENGTH
    export CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_DRIVER
    export CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_URL
    export CROMWELL_BUILD_CENTAUR_PRIOR_SLICK_PROFILE
    export CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT
    export CROMWELL_BUILD_CENTAUR_RESOURCES
    export CROMWELL_BUILD_CENTAUR_SLICK_PROFILE
    export CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS
    export CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY
    export CROMWELL_BUILD_CENTAUR_TYPE
    export CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE
    export CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION
    export CROMWELL_BUILD_CENTAUR_TYPE_STANDARD
    export CROMWELL_BUILD_DOCKER_TAG
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

cromwell::private::install_docker_compose() {
    # Install or upgrade docker-compose so that we get the correct exit codes
    # https://docs.docker.com/compose/release-notes/#1230
    # https://docs.docker.com/compose/install/
    curl \
        --location --fail --silent --show-error \
        "https://github.com/docker/compose/releases/download/1.28.5/docker-compose-$(uname -s)-$(uname -m)" \
        > docker-compose
    sudo mv docker-compose /usr/local/bin
    sudo chmod +x /usr/local/bin/docker-compose
}

cromwell::private::pip_install() {
    local pip_package
    pip_package="${1:?pip_install called without a package}"; shift

    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        sudo -H "$(command -v pip3)" install "${pip_package}" "$@"
    elif [[ "${CROMWELL_BUILD_IS_VIRTUAL_ENV}" == "true" ]]; then
        pip3 install "${pip_package}" "$@"
    else
        pip3 install "${pip_package}" --user "$@"
    fi
}

cromwell::private::upgrade_pip() {
    sudo apt-get install -y python3-pip
    cromwell::private::pip_install pip --upgrade
    cromwell::private::pip_install requests[security] --ignore-installed
}

cromwell::private::install_wait_for_it() {
    curl -s "${CROMWELL_BUILD_WAIT_FOR_IT_URL}" > "$CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT"
    chmod +x "$CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT"
}

cromwell::private::install_vault() {
    curl \
        --location --fail --silent --show-error \
        --output "${CROMWELL_BUILD_VAULT_ZIP}" \
        "https://releases.hashicorp.com/vault/1.6.3/vault_1.6.3_${CROMWELL_BUILD_OS}_amd64.zip"
    unzip "${CROMWELL_BUILD_VAULT_ZIP}" -d "$(dirname "${CROMWELL_BUILD_VAULT_EXECUTABLE}")"
}

cromwell::private::install_git_secrets() {
    # Only install git-secrets on CI. Users should have already installed the executable.
    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        git clone https://github.com/awslabs/git-secrets.git "${CROMWELL_BUILD_GIT_SECRETS_DIRECTORY}"
        pushd "${CROMWELL_BUILD_GIT_SECRETS_DIRECTORY}" > /dev/null
        git checkout "${CROMWELL_BUILD_GIT_SECRETS_COMMIT}"
        export PATH="${PATH}:${PWD}"
        popd > /dev/null
    fi
}

cromwell::private::install_minnie_kenny() {
    # Only install minnie-kenny on CI. Users should have already run the script themselves.
    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        pushd "${CROMWELL_BUILD_ROOT_DIRECTORY}" > /dev/null
        ./minnie-kenny.sh --force
        popd > /dev/null
    fi
}

cromwell::private::start_docker() {
    local docker_image
    local docker_name
    local docker_cid_file
    docker_image="${1:?start_docker called without a docker image}"; shift
    docker_name="$(echo "${docker_image}" | tr "/" "_" | tr ":" "-")_$$"
    docker_cid_file="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${docker_name}.cid"

    docker run --name="${docker_name}" --cidfile="${docker_cid_file}" --detach "$@" "${docker_image}"
    docker logs --follow "${docker_name}" 2>&1 | sed "s/^/$(tput setaf 5)${docker_name}$(tput sgr0) /" &

    cromwell::private::add_exit_function docker rm --force --volumes "$(cat "${docker_cid_file}")"
    cromwell::private::add_exit_function rm "${docker_cid_file}"
}

cromwell::private::start_docker_mysql() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::start_docker_mysql "$@"

    else
        local docker_tag
        local docker_port
        docker_tag="${1?start_docker_mysql called without a docker_tag}"
        docker_port="${2?start_docker_mysql called without a docker_port}"
        shift 2
        cromwell::private::start_docker \
            mysql:"${docker_tag}" \
            --publish "${docker_port}":3306 \
            --env MYSQL_ROOT_PASSWORD=private \
            --env MYSQL_USER="${CROMWELL_BUILD_DATABASE_USERNAME}" \
            --env MYSQL_PASSWORD="${CROMWELL_BUILD_DATABASE_PASSWORD}" \
            --env MYSQL_DATABASE="${CROMWELL_BUILD_DATABASE_SCHEMA}" \
            --volume "${CROMWELL_BUILD_DOCKER_DIRECTORY}"/mysql-conf.d:/etc/mysql/conf.d \

    fi
}

cromwell::private::start_docker_mariadb() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::start_docker_mariadb "$@"

    else
        local docker_tag
        local docker_port
        docker_tag="${1?start_docker_mariadb called without a docker_tag}"
        docker_port="${2?start_docker_mariadb called without a docker_port}"
        shift 2
        cromwell::private::start_docker \
            mariadb:"${docker_tag}" \
            --publish "${docker_port}":3306 \
            --env MYSQL_ROOT_PASSWORD=private \
            --env MYSQL_USER="${CROMWELL_BUILD_DATABASE_USERNAME}" \
            --env MYSQL_PASSWORD="${CROMWELL_BUILD_DATABASE_PASSWORD}" \
            --env MYSQL_DATABASE="${CROMWELL_BUILD_DATABASE_SCHEMA}" \
            --volume "${CROMWELL_BUILD_DOCKER_DIRECTORY}"/mariadb-conf.d:/etc/mysql/conf.d \

    fi
}

cromwell::private::start_docker_postgresql() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::start_docker_postgresql "$@"

    else
        local docker_tag
        local docker_port
        docker_tag="${1?start_docker_postgresql called without a docker_tag}"
        docker_port="${2?start_docker_postgresql called without a docker_port}"
        shift 2
        cromwell::private::start_docker \
            postgres:"${docker_tag}" \
            --publish "${docker_port}":5432 \
            --env POSTGRES_USER="${CROMWELL_BUILD_DATABASE_USERNAME}" \
            --env POSTGRES_PASSWORD="${CROMWELL_BUILD_DATABASE_PASSWORD}" \
            --env POSTGRES_DB="${CROMWELL_BUILD_DATABASE_SCHEMA}" \
            --volume "${CROMWELL_BUILD_DOCKER_DIRECTORY}"/postgresql-initdb.d:/docker-entrypoint-initdb.d \

    fi
}

cromwell::private::start_docker_databases() {
    if [[ -n "${CROMWELL_BUILD_MYSQL_DOCKER_TAG:+set}" ]]; then
        cromwell::private::start_docker_mysql \
            "${CROMWELL_BUILD_MYSQL_DOCKER_TAG}" "${CROMWELL_BUILD_MYSQL_PORT}"
    fi
    if [[ -n "${CROMWELL_BUILD_MARIADB_DOCKER_TAG:+set}" ]]; then
        cromwell::private::start_docker_mariadb \
            "${CROMWELL_BUILD_MARIADB_DOCKER_TAG}" "${CROMWELL_BUILD_MARIADB_PORT}"
    fi
    if [[ -n "${CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG:+set}" ]]; then
        cromwell::private::start_docker_postgresql \
            "${CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG}" "${CROMWELL_BUILD_POSTGRESQL_PORT}"
    fi
    if [[ -n "${CROMWELL_BUILD_MYSQL_LATEST_TAG:+set}" ]]; then
        cromwell::private::start_docker_mysql \
            "${CROMWELL_BUILD_MYSQL_LATEST_TAG}" "${CROMWELL_BUILD_MYSQL_LATEST_PORT}"
    fi
    if [[ -n "${CROMWELL_BUILD_MARIADB_LATEST_TAG:+set}" ]]; then
        cromwell::private::start_docker_mariadb \
            "${CROMWELL_BUILD_MARIADB_LATEST_TAG}" "${CROMWELL_BUILD_MARIADB_LATEST_PORT}"
    fi
    if [[ -n "${CROMWELL_BUILD_POSTGRESQL_LATEST_TAG:+set}" ]]; then
        cromwell::private::start_docker_postgresql \
            "${CROMWELL_BUILD_POSTGRESQL_LATEST_TAG}" "${CROMWELL_BUILD_POSTGRESQL_LATEST_PORT}"
    fi
}

cromwell::private::vault_run() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::vault_run "$@"
    else
        # Run a vault executable that is NOT hosted inside of a docker.io image.
        # For those committers with vault access this avoids pull rate limits reported in BT-143.
        VAULT_ADDR=https://clotho.broadinstitute.org:8200 "${CROMWELL_BUILD_VAULT_EXECUTABLE}" "$@"
    fi
}

cromwell::private::login_vault() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::login_vault
    else
        local vault_token

        # shellcheck disable=SC2153
        if [[ -n "${VAULT_ROLE_ID:+set}" ]] && [[ -n "${VAULT_SECRET_ID:+set}" ]]; then
            vault_token="$(
                cromwell::private::vault_run \
                    write -field=token \
                    auth/approle/login role_id="${VAULT_ROLE_ID}" secret_id="${VAULT_SECRET_ID}"
            )"
        else
            vault_token="${VAULT_TOKEN:-}"
        fi

        if [[ -n "${vault_token}" ]]; then
            # Don't fail here if vault login fails
            # shellcheck disable=SC2015
            cromwell::private::vault_run \
                login "${vault_token}" < /dev/null > /dev/null \
                && echo vault login success \
                || true
        fi
    fi
}

cromwell::private::login_docker() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::login_docker
    else
        local docker_username
        local docker_password

        # Do not fail if docker login fails. We'll try to pull images anonymously.
        docker_username="$(
            cromwell::private::vault_run read -field=username secret/dsde/cromwell/common/cromwell-dockerhub || true
        )"
        docker_password="$(
            cromwell::private::vault_run read -field=password secret/dsde/cromwell/common/cromwell-dockerhub || true
        )"
        docker login --username "${docker_username}" --password-stdin <<< "${docker_password}" || true
    fi
}

cromwell::private::render_secure_resources() {
    # Avoid docker output to sbt's stderr by pulling the image here
    docker pull broadinstitute/dsde-toolbox:dev | cat
    # Copy the CI resources, then render the secure resources using Vault
    sbt -Dsbt.supershell=false --warn renderCiResources \
    || if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
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

cromwell::private::make_build_directories() {
    mkdir -p "${CROMWELL_BUILD_LOG_DIRECTORY}"
    mkdir -p "${CROMWELL_BUILD_RESOURCES_DIRECTORY}"
}

cromwell::private::find_cromwell_jar() {
    CROMWELL_BUILD_CROMWELL_JAR="$( \
        find "${CROMWELL_BUILD_ROOT_DIRECTORY}/server/target/scala-2.13" -name "cromwell-*.jar" -print0 \
        | xargs -0 ls -1 -t \
        | head -n 1 \
        2> /dev/null \
        || true)"
    export CROMWELL_BUILD_CROMWELL_JAR
}

cromwell::private::exists_cromwell_jar() {
    test -s "${CROMWELL_BUILD_CROMWELL_JAR}"
}

cromwell::private::assemble_jars() {
    # CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND allows for an override of the default `assembly` command for assembly.
    # This can be useful to reduce time and memory that might otherwise be spent assembling unused subprojects.
    # shellcheck disable=SC2086
    sbt \
        -Dsbt.supershell=false \
        'set ThisBuild / assembly / logLevel := Level.Error' \
        --warn \
        ${CROMWELL_BUILD_SBT_COVERAGE_COMMAND} \
        --error \
        ${CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND}
}

cromwell::private::setup_prior_version_resources() {
    if [[ "${CROMWELL_BUILD_REQUIRES_PRIOR_VERSION}" == "true" ]]; then
        local prior_config
        local prior_jar
        prior_config="${CROMWELL_BUILD_CROMWELL_CONFIG/%_application.conf/_${CROMWELL_BUILD_PRIOR_VERSION_NUMBER}_application.conf}"
        prior_jar="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_${CROMWELL_BUILD_PRIOR_VERSION_NUMBER}.jar"

        if [[ -f "${prior_config}" ]]; then
            CROMWELL_BUILD_PRE_RESTART_CROMWELL_CONFIG="${prior_config}"
        else
            CROMWELL_BUILD_PRE_RESTART_CROMWELL_CONFIG="${CROMWELL_BUILD_CROMWELL_CONFIG}"
        fi

        CROMWELL_BUILD_PRE_RESTART_DOCKER_TAG="${CROMWELL_BUILD_PRIOR_VERSION_NUMBER}"
        CROMWELL_BUILD_PRE_RESTART_CROMWELL_JAR="${prior_jar}"

        # Copy the prior versions jar out of the previously published docker image
        docker run \
            --rm \
            --entrypoint= \
            --volume "${CROMWELL_BUILD_RESOURCES_DIRECTORY}:${CROMWELL_BUILD_RESOURCES_DIRECTORY}" \
            broadinstitute/cromwell:"${CROMWELL_BUILD_PRE_RESTART_DOCKER_TAG}" \
            cp /app/cromwell.jar "${CROMWELL_BUILD_PRE_RESTART_CROMWELL_JAR}"
    else
        # In tests that are looking for a prior version, actually just use the current version
        CROMWELL_BUILD_PRE_RESTART_CROMWELL_CONFIG="${CROMWELL_BUILD_CROMWELL_CONFIG}"
        CROMWELL_BUILD_PRE_RESTART_DOCKER_TAG="${CROMWELL_BUILD_DOCKER_TAG}"
        CROMWELL_BUILD_PRE_RESTART_CROMWELL_JAR="${CROMWELL_BUILD_CROMWELL_JAR}"
    fi

    export CROMWELL_BUILD_PRE_RESTART_CROMWELL_CONFIG
    export CROMWELL_BUILD_PRE_RESTART_DOCKER_TAG
    export CROMWELL_BUILD_PRE_RESTART_CROMWELL_JAR
}

cromwell::private::generate_code_coverage() {
    sbt -Dsbt.supershell=false --warn coverageReport
    sbt -Dsbt.supershell=false --warn coverageAggregate
    bash <(curl -s https://codecov.io/bash) > /dev/null || true
}

cromwell::private::start_build_heartbeat() {
    # Sleep one minute between printouts, but don't zombie forever
    for ((i=0; i < "${CROMWELL_BUILD_HEARTBEAT_MINUTES}"; i++)); do
        sleep 60
        # shellcheck disable=SC2059
        printf "${CROMWELL_BUILD_HEARTBEAT_PATTERN}"
    done &
    CROMWELL_BUILD_HEARTBEAT_PID=$!
    cromwell::private::add_exit_function cromwell::private::kill_build_heartbeat
}

cromwell::private::start_cromwell_log_tail() {
    while [[ ! -f "${CROMWELL_BUILD_CROMWELL_LOG}" ]]; do
        sleep 2
    done && tail -n 0 -f "${CROMWELL_BUILD_CROMWELL_LOG}" 2> /dev/null &
    CROMWELL_BUILD_CROMWELL_LOG_TAIL_PID=$!
    cromwell::private::add_exit_function cromwell::private::kill_cromwell_log_tail
}

cromwell::private::start_centaur_log_tail() {
    while [[ ! -f "${CROMWELL_BUILD_CENTAUR_LOG}" ]]; do
        sleep 2
    done && tail -n 0 -f "${CROMWELL_BUILD_CENTAUR_LOG}" 2> /dev/null &
    CROMWELL_BUILD_CENTAUR_LOG_TAIL_PID=$!
    cromwell::private::add_exit_function cromwell::private::kill_centaur_log_tail
}

cromwell::private::cat_centaur_log() {
    echo "CENTAUR LOG"
    cat "${CROMWELL_BUILD_CENTAUR_LOG}"
}

cromwell::private::kill_build_heartbeat() {
    if [[ -n "${CROMWELL_BUILD_HEARTBEAT_PID:+set}" ]]; then
        cromwell::private::kill_tree "${CROMWELL_BUILD_HEARTBEAT_PID}"
    fi
}

cromwell::private::kill_cromwell_log_tail() {
    if [[ -n "${CROMWELL_BUILD_CROMWELL_LOG_TAIL_PID:+set}" ]]; then
        cromwell::private::kill_tree "${CROMWELL_BUILD_CROMWELL_LOG_TAIL_PID}"
    fi
}

cromwell::private::kill_centaur_log_tail() {
    if [[ -n "${CROMWELL_BUILD_CENTAUR_LOG_TAIL_PID:+set}" ]]; then
        cromwell::private::kill_tree ${CROMWELL_BUILD_CENTAUR_LOG_TAIL_PID}
    fi
}

cromwell::private::run_exit_functions() {
    if [[ -f "${CROMWELL_BUILD_EXIT_FUNCTIONS}" ]]; then
        local exit_function
        while read -r exit_function; do
          ${exit_function} || true
        done < "${CROMWELL_BUILD_EXIT_FUNCTIONS}"
        rm "${CROMWELL_BUILD_EXIT_FUNCTIONS}" || true
    fi
}

# Adds the function to the list of functions to run on exit.
# Requires at least one positional parameter, the function to run.
cromwell::private::add_exit_function() {
    if [[ "$#" -eq 0 ]]; then
        echo "Error: add_exit_function called without a function" >&2
        exit 1
    fi
    echo "$@" >> "${CROMWELL_BUILD_EXIT_FUNCTIONS}"
    trap cromwell::private::run_exit_functions TERM EXIT
}

cromwell::private::exec_silent_function() {
    local silent_function
    local xtrace_restore_function
    silent_function="${1:?exec_silent_function called without a function}"; shift

    xtrace_restore_function="$(shopt -po xtrace || true)"
    shopt -uo xtrace
    ${silent_function} "$@"
    ${xtrace_restore_function}
}

cromwell::private::is_xtrace_enabled() {
    # Return 0 if xtrace is disabled (set +x), 1 if xtrace is enabled (set -x).
    shopt -qo xtrace
}

cromwell::private::kill_tree() {
  local pid
  local cpid
  pid="${1:?kill_tree called without a pid}"; shift
  for cpid in $(pgrep -P "${pid}"); do
    cromwell::private::kill_tree "${cpid}"
  done
  kill "${pid}" 2> /dev/null
}

cromwell::build::exec_test_script() {
    cromwell::private::create_build_variables
    if [[ "${CROMWELL_BUILD_RUN_TESTS}" == "false" ]]; then
      echo "Use '[force ci]' in commit message to run tests on 'push'"
      exit 0
    fi
    cromwell::private::exec_test_script
}

cromwell::build::setup_common_environment() {
    cromwell::private::check_debug
    cromwell::private::create_build_variables
    cromwell::private::echo_build_variables
    cromwell::private::make_build_directories
    cromwell::private::install_git_secrets
    cromwell::private::install_minnie_kenny
    cromwell::private::install_wait_for_it
    cromwell::private::create_database_variables

    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_GITHUB}")
            # Try to login to vault, and if successful then use vault creds to login to docker.
            # For those committers with vault access this avoids pull rate limits reported in BT-143.
            cromwell::private::install_vault
            cromwell::private::login_vault
            cromwell::private::login_docker
            #Note: Unlike with other CI providers, we are using Github Actions to install Java and sbt for us.
            #This is automatically handled in the set_up_cromwell Github Action, which can be found in
            #[cromwell root]/.github/set_up_cromwell_aciton.
            cromwell::private::install_docker_compose
            cromwell::private::delete_boto_config
            cromwell::private::delete_sbt_boot
            cromwell::private::upgrade_pip
            cromwell::private::start_docker_databases
            ;;
        *)
            ;;
    esac

    cromwell::private::render_secure_resources
    cromwell::private::start_build_heartbeat
}


cromwell::build::setup_centaur_environment() {
    cromwell::private::create_centaur_variables
    cromwell::private::start_cromwell_log_tail
    cromwell::private::start_centaur_log_tail
    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        cromwell::private::add_exit_function cromwell::private::cat_centaur_log
    fi
}

cromwell::private::find_or_assemble_cromwell_jar() {
    cromwell::private::find_cromwell_jar
    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]] || ! cromwell::private::exists_cromwell_jar; then
        echo "Please wait, building jars…"
        cromwell::private::assemble_jars
    fi
    cromwell::private::find_cromwell_jar
    if ! cromwell::private::exists_cromwell_jar; then
        echo "Error: find_or_assemble_cromwell_jar did not locate a cromwell jar even after assembly" >&2
        exit 1
    fi
}

cromwell::build::assemble_jars() {
    cromwell::private::find_or_assemble_cromwell_jar
    cromwell::private::setup_prior_version_resources
}

cromwell::build::build_docker_image() {
    local executable_name
    local docker_image
    executable_name="${1:?build_docker_image called without a executable_name}"
    docker_image="${2:?build_docker_image called without a docker_image}"
    shift
    shift

    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]] || ! docker image ls --quiet "${docker_image}" | grep .; then
        echo "Please wait, building ${executable_name} into ${docker_image}…"

        sbt \
            --error \
            "set \`${executable_name}\`/docker/imageNames := List(ImageName(\"${docker_image}\"))" \
            "${executable_name}/docker"
    fi
}

cromwell::build::build_cromwell_docker() {
    cromwell::build::build_docker_image server broadinstitute/cromwell:"${CROMWELL_BUILD_DOCKER_TAG}"
}

cromwell:build::run_sbt_test() {
    # CROMWELL_BUILD_SBT_COVERAGE_COMMAND allows enabling or disabling `sbt coverage`.
    # Note: sbt logging level now affects the test logging level: https://github.com/sbt/sbt/issues/4480
    # Globally leaving the sbt log level at info for now.
    # Disabling the supershell to reduce log levels.
    # Splitting the JVMs for compilation then scalatest-with-cromwell to reduce memory pressure.
    # Splitting the JVMs for testing-by-sbt-project to also reduce memory pressure.
    # The list of sbt projects is generated by parsing this `log.info()` output, with log color formatting turned off:
    # https://github.com/sbt/sbt/blob/v1.4.9/main/src/main/scala/sbt/Main.scala#L759-L760
    # For more information on testing and memory see also: https://olegych.github.io/blog/sbt-fork.html

    # shellcheck disable=SC2086
    sbt \
        -Dsbt.supershell=false \
        ${CROMWELL_BUILD_SBT_COVERAGE_COMMAND} \
        Test/compile

    local sbt_tests

    if [[ -n "${CROMWELL_BUILD_SBT_INCLUDE}" ]]; then
        # Test only the projects specified
        sbt_tests=$(
            sbt -Dsbt.log.noformat=true projects |
                grep -F $'[info] \t   ' |
                awk '{print $2}' |
                grep -E "^(${CROMWELL_BUILD_SBT_INCLUDE})$" |
                awk '{printf "%s/test ", $1}' \
                || true
        )
    elif [[ -n "${CROMWELL_BUILD_SBT_EXCLUDE}" ]]; then
        # Test all the projects except a few exclusions
        sbt_tests=$(
            sbt -Dsbt.log.noformat=true projects |
                grep -F $'[info] \t   ' |
                awk '{print $2}' |
                grep -v -E "^(${CROMWELL_BUILD_SBT_EXCLUDE})$" |
                awk '{printf "%s/test ", $1}' \
                || true
        )
    else
        # Test all the projects
        sbt_tests="test"
    fi

    # Ensure we are testing something
    if [[ -z "${sbt_tests}" ]]; then
        echo "Error: Unable to retrieve list of sbt projects." >&2
        echo "CROMWELL_BUILD_SBT_INCLUDE='${CROMWELL_BUILD_SBT_INCLUDE}'" >&2
        echo "CROMWELL_BUILD_SBT_EXCLUDE='${CROMWELL_BUILD_SBT_EXCLUDE}'" >&2
        exit 1
    fi

    echo "Starting sbt ${sbt_tests}"
    # shellcheck disable=SC2086
    sbt \
        -Dsbt.supershell=false \
        -Dakka.test.timefactor=${CROMWELL_BUILD_UNIT_SPAN_SCALE_FACTOR} \
        -Dbackend.providers.Local.config.filesystems.local.localization.0=copy \
        ${CROMWELL_BUILD_SBT_COVERAGE_COMMAND} \
        ${sbt_tests}
}

cromwell::build::run_centaur() {
    local -a additional_args
    additional_args=()
    if [[ -n "${CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS-}" ]]; then
        # Allow splitting on space to simulate an exported array
        # https://stackoverflow.com/questions/5564418/exporting-an-array-in-bash-script#answer-5564589
        # shellcheck disable=SC2206
        additional_args=(${CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS})
    fi
    if [[ "${CROMWELL_BUILD_GENERATE_COVERAGE}" == "true" ]]; then
        additional_args+=("-g")
    fi
    # Handle empty arrays in older versions of bash
    # https://stackoverflow.com/questions/7577052/bash-empty-array-expansion-with-set-u#answer-7577209
    "${CROMWELL_BUILD_ROOT_DIRECTORY}/centaur/test_cromwell.sh" \
        -n "${CROMWELL_BUILD_CENTAUR_CONFIG}" \
        -l "${CROMWELL_BUILD_LOG_DIRECTORY}" \
        ${additional_args[@]+"${additional_args[@]}"} \
        "$@"
}

cromwell::build::generate_code_coverage() {
    if [[ "${CROMWELL_BUILD_GENERATE_COVERAGE}" == "true" ]]; then
        cromwell::private::generate_code_coverage
    fi
}

cromwell::build::print_workflow_statistics() {
    echo "Total workflows"
    mysql --host=127.0.0.1 --user=cromwell --password=test cromwell_test -e \
        "SELECT COUNT(*) as total_workflows_run FROM WORKFLOW_METADATA_SUMMARY_ENTRY;"

    echo "Late starters"
    mysql --host=127.0.0.1 --user=cromwell --password=test cromwell_test -e \
        "SELECT WORKFlOW_NAME as name,
            TIMESTAMPDIFF(MINUTE, START_TIMESTAMP, END_TIMESTAMP) as runtime_minutes,
            START_TIMESTAMP as START,
            END_TIMESTAMP as end
        FROM WORKFLOW_METADATA_SUMMARY_ENTRY
            WHERE PARENT_WORKFLOW_EXECUTION_UUID IS NULL # exclude subworkflows
            ORDER BY START_TIMESTAMP DESC
            LIMIT 20;"

    echo "Late finishers"
    mysql --host=127.0.0.1 --user=cromwell --password=test cromwell_test -e \
        "SELECT WORKFlOW_NAME as name,
            TIMESTAMPDIFF(MINUTE, START_TIMESTAMP, END_TIMESTAMP) as runtime_minutes,
            START_TIMESTAMP as start,
            END_TIMESTAMP as END
        FROM WORKFLOW_METADATA_SUMMARY_ENTRY
            WHERE PARENT_WORKFLOW_EXECUTION_UUID IS NULL
            ORDER BY END_TIMESTAMP DESC
            LIMIT 20;"

    echo "Long duration"
    mysql --host=127.0.0.1 --user=cromwell --password=test cromwell_test -e \
        "SELECT WORKFlOW_NAME as name,
            TIMESTAMPDIFF(MINUTE, START_TIMESTAMP, END_TIMESTAMP) as RUNTIME_MINUTES,
            START_TIMESTAMP as start,
            END_TIMESTAMP as end
        FROM WORKFLOW_METADATA_SUMMARY_ENTRY
            WHERE PARENT_WORKFLOW_EXECUTION_UUID IS NULL
            ORDER BY RUNTIME_MINUTES DESC
            LIMIT 20;"
}

cromwell::build::exec_silent_function() {
    local silent_function
    silent_function="${1:?exec_silent_function called without a function}"; shift
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function "${silent_function}" "$@"
    else
        ${silent_function} "$@"
    fi
}

cromwell::build::pip_install() {
    cromwell::private::pip_install "$@"
}

cromwell::build::add_exit_function() {
    cromwell::private::add_exit_function "$1"
}

cromwell::build::delete_docker_images() {
    local docker_delete_function
    local docker_image_file
    docker_delete_function="${1:?delete_images called without a docker_delete_function}"
    docker_image_file="${2:?delete_images called without a docker_image_file}"
    shift
    shift

    if [[ -f "${docker_image_file}" ]]; then
        local docker_image
        while read -r docker_image; do
          ${docker_delete_function} "${docker_image}" || true
        done < "${docker_image_file}"
        rm "${docker_image_file}" || true
    fi
}

cromwell::build::kill_tree() {
    cromwell::private::kill_tree "$1"
}
