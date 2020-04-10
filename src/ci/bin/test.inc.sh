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

cromwell::private::check_debug() {
    # shellcheck disable=SC2154
    if [[ -n "${crmdbg:+set}" ]]; then
        set -o xtrace
    fi

    # shellcheck disable=SC2154
    if [[ -n "${crmcit:+set}" ]]; then
        CROMWELL_BUILD_CENTAUR_TYPE="integration"
    fi
}

# Exports environment variables used for scripts.
cromwell::private::create_build_variables() {
    CROMWELL_BUILD_PROVIDER_TRAVIS="travis"
    CROMWELL_BUILD_PROVIDER_JENKINS="jenkins"
    CROMWELL_BUILD_PROVIDER_UNKNOWN="unknown"

    if [[ "${TRAVIS-false}" == "true" ]]; then
        CROMWELL_BUILD_PROVIDER="${CROMWELL_BUILD_PROVIDER_TRAVIS}"
    elif [[ "${JENKINS-false}" == "true" ]]; then
        CROMWELL_BUILD_PROVIDER="${CROMWELL_BUILD_PROVIDER_JENKINS}"
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
        CROMWELL_BUILD_PRIOR_VERSION_NUMBER=$((${CROMWELL_BUILD_CURRENT_VERSION_NUMBER} - 1))
    fi

    local git_commit_message
    git_commit_message="$(git log --format=%B --max-count=1 HEAD 2>/dev/null || true)"
    if [[ "${git_commit_message}" == *"[force ci]"* ]]; then
        CROMWELL_BUILD_FORCE_TESTS=true
    else
        CROMWELL_BUILD_FORCE_TESTS=false
    fi

    local git_revision
    if git_revision="$(git rev-parse --short=7 HEAD 2>/dev/null)"; then
        CROMWELL_BUILD_GIT_HASH_SUFFIX="g${git_revision}"
    else
        CROMWELL_BUILD_GIT_HASH_SUFFIX="gUNKNOWN"
    fi

    # Value of the `TRAVIS_BRANCH` variable depends on type of Travis build: if it is pull request build, the value
    # will be the name of the branch targeted by the pull request, and for push builds it will be the name of the
    # branch. So, in case of push builds `git diff` will always return empty result. This is why we only use this short
    # circuiting logic for pull request builds
    if [[ "${TRAVIS_EVENT_TYPE:-unset}" != "pull_request" ]]; then
        CROMWELL_BUILD_ONLY_DOCS_CHANGED=false
    elif git diff --name-only "origin/${TRAVIS_BRANCH}" 2>&1 | egrep -q --invert-match "^mkdocs.yml|^docs/"; then
        CROMWELL_BUILD_ONLY_DOCS_CHANGED=false
    else
        CROMWELL_BUILD_ONLY_DOCS_CHANGED=true
    fi

    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
            CROMWELL_BUILD_IS_CI=true
            CROMWELL_BUILD_IS_SECURE="${TRAVIS_SECURE_ENV_VARS}"

            CROMWELL_BUILD_TYPE="${BUILD_TYPE}"
            CROMWELL_BUILD_BRANCH="${TRAVIS_PULL_REQUEST_BRANCH:-${TRAVIS_BRANCH}}"
            CROMWELL_BUILD_EVENT="${TRAVIS_EVENT_TYPE}"
            CROMWELL_BUILD_TAG="${TRAVIS_TAG}"
            CROMWELL_BUILD_NUMBER="${TRAVIS_JOB_NUMBER}"
            CROMWELL_BUILD_URL="https://travis-ci.com/${TRAVIS_REPO_SLUG}/jobs/${TRAVIS_JOB_ID}"
            CROMWELL_BUILD_GIT_USER_EMAIL="travis@travis-ci.com"
            CROMWELL_BUILD_GIT_USER_NAME="Travis CI"
            CROMWELL_BUILD_HEARTBEAT_PATTERN="…"
            CROMWELL_BUILD_GENERATE_COVERAGE=true

            # For solely documentation updates run only checkPublish. Otherwise always run sbt, even for 'push'.
            # This allows quick sanity checks before starting PRs *and* publishing after merges into develop.
            if [[ "${CROMWELL_BUILD_FORCE_TESTS}" == "true" ]]; then
                CROMWELL_BUILD_RUN_TESTS=true
            elif [[ "${CROMWELL_BUILD_ONLY_DOCS_CHANGED}" == "true" ]] && \
                [[ "${BUILD_TYPE}" != "checkPublish" ]]; then
                CROMWELL_BUILD_RUN_TESTS=false
            elif [[ "${TRAVIS_EVENT_TYPE}" == "push" ]] && \
                [[ "${BUILD_TYPE}" != "sbt" ]]; then
                CROMWELL_BUILD_RUN_TESTS=false
            else
                CROMWELL_BUILD_RUN_TESTS=true
            fi
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
            CROMWELL_BUILD_HEARTBEAT_PATTERN="…\n"
            CROMWELL_BUILD_GENERATE_COVERAGE=false
            CROMWELL_BUILD_RUN_TESTS=true
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
            CROMWELL_BUILD_HEARTBEAT_PATTERN="…"
            CROMWELL_BUILD_GENERATE_COVERAGE=true
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
    backend_type="${backend_type#conformance}"
    backend_type="$(echo "${backend_type}" | sed 's/\([A-Z]\)/_\1/g' | tr '[:upper:]' '[:lower:]' | cut -c 2-)"
    CROMWELL_BUILD_BACKEND_TYPE="${backend_type}"

    if [[ "${CROMWELL_BUILD_TYPE}" == conformance* ]]; then
        CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND="server/assembly centaurCwlRunner/assembly"
    else
        CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND="assembly"
    fi

    case "${CROMWELL_BUILD_TYPE}" in
        centaurPapiUpgradePapiV1*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v1_v2alpha1_upgrade_application.conf"
            ;;
        centaurPapiUpgradeNewWorkflowsPapiV1*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v1_v2alpha1_upgrade_application.conf"
            ;;
        centaurPapiUpgradePapiV2alpha1*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2alpha1_v2beta_upgrade_application.conf"
            ;;
        centaurPapiUpgradeNewWorkflowsPapiV2alpha1*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2alpha1_v2beta_upgrade_application.conf"
            ;;
        centaurHoricromtalPapiV2alpha1*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2alpha1_horicromtal_application.conf"
            ;;
        centaurHoricromtalPapiV2beta*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2beta_horicromtal_application.conf"
            ;;
        centaurHoricromtalEngineUpgrade*)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/papi_v2alpha1_horicromtal_application.conf"
            ;;
        *)
            CROMWELL_BUILD_CROMWELL_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${CROMWELL_BUILD_BACKEND_TYPE}_application.conf"
            ;;
    esac

    CROMWELL_BUILD_OPTIONAL_SECURE="${CROMWELL_BUILD_OPTIONAL_SECURE-false}"
    CROMWELL_BUILD_REQUIRES_SECURE="${CROMWELL_BUILD_REQUIRES_SECURE-false}"
    CROMWELL_BUILD_REQUIRES_PRIOR_VERSION="${CROMWELL_BUILD_REQUIRES_PRIOR_VERSION-false}"
    CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND="${CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND-assembly}"
    CROMWELL_BUILD_CROMWELL_DOCKER_TAG="${CROMWELL_BUILD_CROMWELL_DOCKER_TAG-just-testing-horicromtal}"
    CROMWELL_BUILD_CROMWELL_DOCKER_COMPOSE="${CROMWELL_BUILD_DOCKER_DIRECTORY}/docker-compose-horicromtal.yml"
    VAULT_TOKEN="${VAULT_TOKEN-vault token is not set as an environment variable}"

    local hours_to_minutes
    hours_to_minutes=60
    CROMWELL_BUILD_HEARTBEAT_MINUTES=$((20 * hours_to_minutes))

    export CROMWELL_BUILD_UNIT_TEST_EXCLUDE_TAGS
    export CROMWELL_BUILD_BACKEND_TYPE
    export CROMWELL_BUILD_BRANCH
    export CROMWELL_BUILD_CROMWELL_CONFIG
    export CROMWELL_BUILD_CROMWELL_LOG
    export CROMWELL_BUILD_CURRENT_VERSION_NUMBER
    export CROMWELL_BUILD_DOCKER_DIRECTORY
    export CROMWELL_BUILD_EVENT
    export CROMWELL_BUILD_EXIT_FUNCTIONS
    export CROMWELL_BUILD_FORCE_TESTS
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
    export CROMWELL_BUILD_IS_SECURE
    export CROMWELL_BUILD_IS_VIRTUAL_ENV
    export CROMWELL_BUILD_LOG_DIRECTORY
    export CROMWELL_BUILD_NUMBER
    export CROMWELL_BUILD_OPTIONAL_SECURE
    export CROMWELL_BUILD_OS
    export CROMWELL_BUILD_OS_DARWIN
    export CROMWELL_BUILD_OS_LINUX
    export CROMWELL_BUILD_PRIOR_VERSION_NUMBER
    export CROMWELL_BUILD_PROVIDER
    export CROMWELL_BUILD_PROVIDER_JENKINS
    export CROMWELL_BUILD_PROVIDER_TRAVIS
    export CROMWELL_BUILD_PROVIDER_UNKNOWN
    export CROMWELL_BUILD_REQUIRES_SECURE
    export CROMWELL_BUILD_REQUIRES_PRIOR_VERSION
    export CROMWELL_BUILD_RESOURCES_DIRECTORY
    export CROMWELL_BUILD_RESOURCES_SOURCES
    export CROMWELL_BUILD_ROOT_DIRECTORY
    export CROMWELL_BUILD_RUN_TESTS
    export CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND
    export CROMWELL_BUILD_SCRIPTS_DIRECTORY
    export CROMWELL_BUILD_TAG
    export CROMWELL_BUILD_CROMWELL_DOCKER_COMPOSE
    export CROMWELL_BUILD_CROMWELL_DOCKER_TAG
    export CROMWELL_BUILD_TYPE
    export CROMWELL_BUILD_URL
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
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
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
        "${CROMWELL_BUILD_PROVIDER_JENKINS}")
            # NOTE: Jenkins uses src/ci/docker-compose/docker-compose.yml.
            # We don't define a docker tag because the docker-compose has already spun up the database containers by the
            # time this script is run. Other variables here must match the database service names and settings the yaml.
            CROMWELL_BUILD_MARIADB_HOSTNAME="mariadb-db"
            CROMWELL_BUILD_MARIADB_PORT="3306"
            CROMWELL_BUILD_MARIADB_DOCKER_TAG=""
            CROMWELL_BUILD_MARIADB_LATEST_HOSTNAME="mariadb-db-latest"
            CROMWELL_BUILD_MARIADB_LATEST_PORT="3306"
            CROMWELL_BUILD_MARIADB_LATEST_TAG=""
            CROMWELL_BUILD_MYSQL_HOSTNAME="mysql-db"
            CROMWELL_BUILD_MYSQL_PORT="3306"
            CROMWELL_BUILD_MYSQL_DOCKER_TAG=""
            CROMWELL_BUILD_MYSQL_LATEST_HOSTNAME="mysql-db-latest"
            CROMWELL_BUILD_MYSQL_LATEST_PORT="3306"
            CROMWELL_BUILD_MYSQL_LATEST_TAG=""
            CROMWELL_BUILD_POSTGRESQL_HOSTNAME="postgresql-db"
            CROMWELL_BUILD_POSTGRESQL_PORT="5432"
            CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG=""
            CROMWELL_BUILD_POSTGRESQL_LATEST_HOSTNAME="postgresql-db-latest"
            CROMWELL_BUILD_POSTGRESQL_LATEST_PORT="3306"
            CROMWELL_BUILD_POSTGRESQL_LATEST_TAG=""
            ;;
        *)
            CROMWELL_BUILD_MARIADB_HOSTNAME="${CROMWELL_BUILD_MARIADB_HOSTNAME-localhost}"
            CROMWELL_BUILD_MARIADB_PORT="${CROMWELL_BUILD_MARIADB_PORT-13306}"
            CROMWELL_BUILD_MARIADB_DOCKER_TAG=""
            CROMWELL_BUILD_MARIADB_LATEST_HOSTNAME="${CROMWELL_BUILD_MARIADB_LATEST_HOSTNAME-localhost}"
            CROMWELL_BUILD_MARIADB_LATEST_PORT="${CROMWELL_BUILD_MARIADB_LATEST_PORT-13306}"
            CROMWELL_BUILD_MARIADB_LATEST_TAG=""
            CROMWELL_BUILD_MYSQL_HOSTNAME="${CROMWELL_BUILD_MYSQL_HOSTNAME-localhost}"
            CROMWELL_BUILD_MYSQL_PORT="${CROMWELL_BUILD_MYSQL_PORT-3306}"
            CROMWELL_BUILD_MYSQL_DOCKER_TAG=""
            CROMWELL_BUILD_MYSQL_LATEST_HOSTNAME="${CROMWELL_BUILD_MYSQL_LATEST_HOSTNAME-localhost}"
            CROMWELL_BUILD_MYSQL_LATEST_PORT="${CROMWELL_BUILD_MYSQL_LATEST_PORT-13306}"
            CROMWELL_BUILD_MYSQL_LATEST_TAG=""
            CROMWELL_BUILD_POSTGRESQL_HOSTNAME="${CROMWELL_BUILD_POSTGRESQL_HOSTNAME-localhost}"
            CROMWELL_BUILD_POSTGRESQL_PORT="${CROMWELL_BUILD_POSTGRESQL_PORT-5432}"
            CROMWELL_BUILD_POSTGRESQL_DOCKER_TAG=""
            CROMWELL_BUILD_POSTGRESQL_LATEST_HOSTNAME="${CROMWELL_BUILD_POSTGRESQL_LATEST_HOSTNAME-localhost}"
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

            # Determine horicromtal Centaur config:
            if test "${CROMWELL_BUILD_BACKEND_TYPE}" = "papi_v2alpha1" || test "${CROMWELL_BUILD_BACKEND_TYPE}" = "papi_v2beta"
            then
              CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application_papi_v2_horicromtal.conf"
            else
              CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application_horicromtal.conf"
            fi
            echo "*** Using centaur config '${CROMWELL_BUILD_CENTAUR_CONFIG}' for backend type '${CROMWELL_BUILD_BACKEND_TYPE}'"
            ;;
        "${CROMWELL_BUILD_CENTAUR_TYPE_HORICROMTAL_ENGINE_UPGRADE}")
            # Use the engine upgrade test cases despite the horicromtal Centaur build type.
            CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY="${CROMWELL_BUILD_CENTAUR_RESOURCES}/engineUpgradeTestCases"
            # Special horicromtal engine upgrade Centaur config with horicromtal assertions turned off.
            CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application_horicromtal_no_assert.conf"
            ;;
        *)
            CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY="${CROMWELL_BUILD_CENTAUR_RESOURCES}/${CROMWELL_BUILD_CENTAUR_TYPE}TestCases"
            if test "${CROMWELL_BUILD_BACKEND_TYPE}" = "papi_v2alpha1" || test "${CROMWELL_BUILD_BACKEND_TYPE}" = "papi_v2beta"
            then
              CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application_papi_v2.conf"
            else
              CROMWELL_BUILD_CENTAUR_CONFIG="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/centaur_application.conf"
            fi
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
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")

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
        "${CROMWELL_BUILD_PROVIDER_JENKINS}")
            CROMWELL_BUILD_CENTAUR_SLICK_PROFILE="slick.jdbc.MySQLProfile$"
            CROMWELL_BUILD_CENTAUR_JDBC_DRIVER="com.mysql.cj.jdbc.Driver"
            CROMWELL_BUILD_CENTAUR_JDBC_URL="${mysql_jdbc_url}"
            CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS="${CENTAUR_TEST_ADDITIONAL_PARAMETERS-}"
            ;;
        *)
            CROMWELL_BUILD_CENTAUR_SLICK_PROFILE="${CROMWELL_BUILD_CENTAUR_SLICK_PROFILE-slick.jdbc.MySQLProfile\$}"
            CROMWELL_BUILD_CENTAUR_JDBC_DRIVER="${CROMWELL_BUILD_CENTAUR_JDBC_DRIVER-com.mysql.cj.jdbc.Driver}"
            CROMWELL_BUILD_CENTAUR_JDBC_URL="${CROMWELL_BUILD_CENTAUR_JDBC_URL-${mysql_jdbc_url}}"
            CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS=
            ;;
    esac

    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        CROMWELL_BUILD_CENTAUR_DOCKER_TAG="${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_NUMBER}"
    else
        CROMWELL_BUILD_CENTAUR_DOCKER_TAG="${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_TYPE}-${CROMWELL_BUILD_GIT_HASH_SUFFIX}"
    fi

    # Trim and replace invalid characters in the docker tag
    # https://docs.docker.com/engine/reference/commandline/tag/#extended-description
    CROMWELL_BUILD_CENTAUR_DOCKER_TAG="${CROMWELL_BUILD_CENTAUR_DOCKER_TAG:0:128}"
    CROMWELL_BUILD_CENTAUR_DOCKER_TAG="${CROMWELL_BUILD_CENTAUR_DOCKER_TAG//[^a-zA-Z0-9.-]/_}"

    case "${CROMWELL_BUILD_CENTAUR_TYPE}" in
        "${CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION}")
            CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=512000
            ;;
        *)
            CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT=128000
            ;;
    esac

    # When upgrading to the MariaDB driver, start with MySQL then switch to MariaDB.
    if [[ "${CROMWELL_BUILD_PROVIDER}" == "${CROMWELL_BUILD_PROVIDER_TRAVIS}" ]] && \
        [[ -n "${CROMWELL_BUILD_MARIADB_DOCKER_TAG:+set}" ]]; then

        CROMWELL_BUILD_CENTAUR_PRIOR_SLICK_PROFILE="slick.jdbc.MySQLProfile$"
        CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_DRIVER="com.mysql.cj.jdbc.Driver"
        CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_URL="jdbc:mysql://${CROMWELL_BUILD_MARIADB_HOSTNAME}:${CROMWELL_BUILD_MARIADB_PORT}/${CROMWELL_BUILD_DATABASE_SCHEMA}?allowPublicKeyRetrieval=true&useSSL=false&rewriteBatchedStatements=true&serverTimezone=UTC&useInformationSchema=true"
    else

        CROMWELL_BUILD_CENTAUR_PRIOR_SLICK_PROFILE="${CROMWELL_BUILD_CENTAUR_PRIOR_SLICK_PROFILE-${CROMWELL_BUILD_CENTAUR_SLICK_PROFILE}}"
        CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_DRIVER="${CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_DRIVER-${CROMWELL_BUILD_CENTAUR_JDBC_DRIVER}}"
        CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_URL="${CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_URL-${CROMWELL_BUILD_CENTAUR_JDBC_URL}}"
    fi

    CROMWELL_BUILD_CENTAUR_256_BITS_KEY="$(dd bs=1 count=32 if=/dev/urandom 2> /dev/null | base64 | tr -d '\n')"

    export CROMWELL_BUILD_CENTAUR_256_BITS_KEY
    export CROMWELL_BUILD_CENTAUR_CONFIG
    export CROMWELL_BUILD_CENTAUR_DOCKER_TAG
    export CROMWELL_BUILD_CENTAUR_JDBC_DRIVER
    export CROMWELL_BUILD_CENTAUR_JDBC_PASSWORD
    export CROMWELL_BUILD_CENTAUR_JDBC_URL
    export CROMWELL_BUILD_CENTAUR_JDBC_USERNAME
    export CROMWELL_BUILD_CENTAUR_LOG
    export CROMWELL_BUILD_CENTAUR_TEST_ADDITIONAL_PARAMETERS
    export CROMWELL_BUILD_CENTAUR_TEST_DIRECTORY
    export CROMWELL_BUILD_CENTAUR_READ_LINES_LIMIT
    export CROMWELL_BUILD_CENTAUR_RESOURCES
    export CROMWELL_BUILD_CENTAUR_SLICK_PROFILE
    export CROMWELL_BUILD_CENTAUR_TYPE
    export CROMWELL_BUILD_CENTAUR_TYPE_STANDARD
    export CROMWELL_BUILD_CENTAUR_TYPE_INTEGRATION
    export CROMWELL_BUILD_CENTAUR_TYPE_ENGINE_UPGRADE
    export CROMWELL_BUILD_CENTAUR_PRIOR_SLICK_PROFILE
    export CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_DRIVER
    export CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_USERNAME
    export CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_PASSWORD
    export CROMWELL_BUILD_CENTAUR_PRIOR_JDBC_URL
}

cromwell::private::create_conformance_variables() {
    CROMWELL_BUILD_CWL_RUNNER_MODE="${CROMWELL_BUILD_BACKEND_TYPE}"
    CROMWELL_BUILD_CWL_TOOL_VERSION="1.0.20190228155703"
    CROMWELL_BUILD_CWL_TEST_VERSION="1.0.20190228134645"
    CROMWELL_BUILD_CWL_TEST_COMMIT="1f501e38ff692a408e16b246ac7d64d32f0822c2" # use known git hash to avoid changes
    CROMWELL_BUILD_CWL_TEST_RUNNER="${CROMWELL_BUILD_ROOT_DIRECTORY}/centaurCwlRunner/src/bin/centaur-cwl-runner.bash"
    CROMWELL_BUILD_CWL_TEST_DIRECTORY="${CROMWELL_BUILD_ROOT_DIRECTORY}/common-workflow-language"
    CROMWELL_BUILD_CWL_TEST_RESOURCES="${CROMWELL_BUILD_CWL_TEST_DIRECTORY}/v1.0/v1.0"
    CROMWELL_BUILD_CWL_TEST_WDL="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cwl_conformance_test.wdl"
    CROMWELL_BUILD_CWL_TEST_INPUTS="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cwl_conformance_test.inputs.json"
    CROMWELL_BUILD_CWL_TEST_OUTPUT="${CROMWELL_BUILD_LOG_DIRECTORY}/cwl_conformance_test.out.txt"
    CROMWELL_BUILD_CWL_TEST_PARALLELISM=10 # Set too high will cause false negatives due to cromwell server timeouts.

    export CROMWELL_BUILD_CWL_RUNNER_MODE
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

cromwell::private::verify_secure_build() {
    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
            if [[ "${CROMWELL_BUILD_IS_SECURE}" != "true" ]] && \
                [[ "${CROMWELL_BUILD_REQUIRES_SECURE}" == "true" ]]; then
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

cromwell::private::exec_test_script() {
    local upper_build_type
    upper_build_type="$(tr '[:lower:]' '[:upper:]' <<< "${CROMWELL_BUILD_TYPE:0:1}")${CROMWELL_BUILD_TYPE:1}"
    exec "${CROMWELL_BUILD_SCRIPTS_DIRECTORY}/test${upper_build_type}.sh"
}

cromwell::private::stop_travis_defaults() {
  # https://stackoverflow.com/questions/27382295/how-to-stop-services-on-travis-ci-running-by-default#answer-27410479
  sudo /etc/init.d/mysql stop
  sudo /etc/init.d/postgresql stop
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

cromwell::private::pip_install() {
    local pip_package
    pip_package="${1:?pip_install called without a package}"; shift

    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        sudo -H pip install "${pip_package}" "$@"
    elif [[ "${CROMWELL_BUILD_IS_VIRTUAL_ENV}" == "true" ]]; then
        pip install "${pip_package}" "$@"
    else
        pip install "${pip_package}" --user "$@"
    fi
}

cromwell::private::upgrade_pip() {
    cromwell::private::pip_install pip --upgrade
    cromwell::private::pip_install requests[security] --ignore-installed
}

cromwell::private::install_wait_for_it() {
    curl -s "${CROMWELL_BUILD_WAIT_FOR_IT_URL}" > "$CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT"
    chmod +x "$CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT"
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

cromwell::private::pull_common_docker_images() {
    # All tests use ubuntu:latest - make sure it's there before starting the tests
    # because pulling the image during some of the tests would cause them to fail
    # (specifically output_redirection which expects a specific value in stderr)
    docker pull ubuntu
}

cromwell::private::install_cwltest() {
    # TODO: No clue why these are needed for cwltool. If you know please update this comment.
    sudo apt-get install procps || true
    cromwell::private::pip_install cwltool=="${CROMWELL_BUILD_CWL_TOOL_VERSION}" --ignore-installed
    cromwell::private::pip_install cwltest=="${CROMWELL_BUILD_CWL_TEST_VERSION}"
}

cromwell::private::checkout_pinned_cwl() {
    if [[ ! -d "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}" ]]; then
        git clone \
            https://github.com/common-workflow-language/common-workflow-language.git \
            "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}"
        (
            pushd "${CROMWELL_BUILD_CWL_TEST_DIRECTORY}" > /dev/null
            git checkout "${CROMWELL_BUILD_CWL_TEST_COMMIT}"
            popd > /dev/null
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
        if [[ -f "${dockerhub_auth_include}" ]]; then
            # shellcheck source=/dev/null
            source "${dockerhub_auth_include}"
        fi
    fi
}

cromwell::private::vault_login() {
    if cromwell::private::is_xtrace_enabled; then
        cromwell::private::exec_silent_function cromwell::private::vault_login
    elif [[ "${CROMWELL_BUILD_IS_SECURE}" == "true" ]]; then
        case "${CROMWELL_BUILD_PROVIDER}" in
            "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
                # Login to vault to access secrets
                local vault_token
                vault_token="${VAULT_TOKEN}"
                # Don't fail here if vault login fails
                # shellcheck disable=SC2015
                docker run --rm \
                    -v "${CROMWELL_BUILD_HOME_DIRECTORY}:/root:rw" \
                    broadinstitute/dsde-toolbox:dev \
                    vault auth "${vault_token}" < /dev/null > /dev/null && echo vault auth success \
                || true
                ;;
            *)
                ;;
        esac
    fi
}

cromwell::private::render_secure_resources() {
    # Copy the CI resources, then render the secure resources using Vault
    sbt --warn renderCiResources \
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

cromwell::private::copy_all_resources() {
    # Only copy the CI resources. Secure resources are not rendered.
    sbt --warn copyCiResources
}

cromwell::private::setup_secure_resources() {
    if [[ "${CROMWELL_BUILD_REQUIRES_SECURE}" == "true" ]] || [[ "${CROMWELL_BUILD_OPTIONAL_SECURE}" == "true" ]]; then
        case "${CROMWELL_BUILD_PROVIDER}" in
            "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
                cromwell::private::vault_login
                cromwell::private::render_secure_resources
                cromwell::private::docker_login
                ;;
            "${CROMWELL_BUILD_PROVIDER_JENKINS}")
                cromwell::private::copy_all_resources
                ;;
            *)
                cromwell::private::render_secure_resources
                ;;
        esac
    else
        cromwell::private::copy_all_resources
    fi
}

cromwell::private::make_build_directories() {
    if [[ "${CROMWELL_BUILD_PROVIDER}" == "${CROMWELL_BUILD_PROVIDER_JENKINS}" ]]; then
        sudo chmod -R a+w .
    fi
    mkdir -p "${CROMWELL_BUILD_LOG_DIRECTORY}"
    mkdir -p "${CROMWELL_BUILD_RESOURCES_DIRECTORY}"
}

cromwell::private::find_cromwell_jar() {
    CROMWELL_BUILD_CROMWELL_JAR="$( \
        find "${CROMWELL_BUILD_ROOT_DIRECTORY}/server/target/scala-2.12" -name "cromwell-*.jar" -print0 \
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
    CROMWELL_SBT_ASSEMBLY_LOG_LEVEL=error sbt --warn coverage ${CROMWELL_BUILD_SBT_ASSEMBLY_COMMAND} -error
}

cromwell::private::setup_prior_version_resources() {
    if [[ "${CROMWELL_BUILD_REQUIRES_PRIOR_VERSION}" == "true" ]]; then
        local prior_config
        local prior_jar
        prior_config="${CROMWELL_BUILD_CROMWELL_CONFIG/%_application.conf/_${CROMWELL_BUILD_PRIOR_VERSION_NUMBER}_application.conf}"
        prior_jar="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_${CROMWELL_BUILD_PRIOR_VERSION_NUMBER}.jar"

        if [[ -f "${prior_config}" ]]; then
            CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_CONFIG="${prior_config}"
        else
            CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_CONFIG="${CROMWELL_BUILD_CROMWELL_CONFIG}"
        fi

        CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_DOCKER_TAG="${CROMWELL_BUILD_PRIOR_VERSION_NUMBER}"
        CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_JAR="${prior_jar}"

        # Copy the prior versions jar out of the previously published docker image
        docker run \
            --rm \
            --entrypoint= \
            --volume "${CROMWELL_BUILD_RESOURCES_DIRECTORY}:${CROMWELL_BUILD_RESOURCES_DIRECTORY}" \
            broadinstitute/cromwell:"${CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_DOCKER_TAG}" \
            cp /app/cromwell.jar "${CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_JAR}"
    else
        # In tests that are looking for a prior version, actually just use the current version
        CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_CONFIG="${CROMWELL_BUILD_CROMWELL_CONFIG}"
        CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_DOCKER_TAG="${CROMWELL_BUILD_CROMWELL_DOCKER_TAG}"
        CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_JAR="${CROMWELL_BUILD_CROMWELL_JAR}"
    fi

    export CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_CONFIG
    export CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_DOCKER_TAG
    export CROMWELL_BUILD_CROMWELL_PRIOR_VERSION_JAR
}

cromwell::private::exists_cromwell_docker() {
    docker image ls --quiet broadinstitute/cromwell:"${CROMWELL_BUILD_CROMWELL_DOCKER_TAG}" | grep .
}

cromwell::private::build_cromwell_docker() {
    CROMWELL_SBT_DOCKER_TAGS="${CROMWELL_BUILD_CROMWELL_DOCKER_TAG}" sbt server/docker -error
}

cromwell::private::generate_code_coverage() {
    sbt --warn coverageReport -warn
    sbt --warn coverageAggregate -warn
    bash <(curl -s https://codecov.io/bash) > /dev/null || true
}

cromwell::private::publish_artifacts_only() {
    CROMWELL_SBT_ASSEMBLY_LOG_LEVEL=warn sbt "$@" publish -warn
}

cromwell::private::publish_artifacts_and_docker() {
    CROMWELL_SBT_ASSEMBLY_LOG_LEVEL=warn sbt "$@" publish dockerBuildAndPush -warn
}

cromwell::private::publish_artifacts_check() {
    sbt --warn verifyArtifactoryCredentialsExist -warn
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
    pushd publish_complete > /dev/null

    git init
    git config core.sshCommand "ssh -i ${github_private_deploy_key} -F /dev/null"
    git config user.email "${CROMWELL_BUILD_GIT_USER_EMAIL}"
    git config user.name "${CROMWELL_BUILD_GIT_USER_NAME}"

    git remote add "${git_publish_remote}" "${git_repo}"
    git checkout -b "${git_publish_branch}"
    git commit --allow-empty -m "${git_publish_message}"
    git push -f "${git_publish_remote}" "${git_publish_branch}"

    popd > /dev/null
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

cromwell::private::cat_conformance_log() {
    echo "CONFORMANCE LOG"
    cat "${CROMWELL_BUILD_CWL_TEST_OUTPUT}"
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

cromwell::private::start_conformance_cromwell() {
    # Start the Cromwell server in the directory containing input files so it can access them via their relative path
    pushd "${CROMWELL_BUILD_CWL_TEST_RESOURCES}" > /dev/null

    # Turn off call caching as hashing doesn't work since it sees local and not GCS paths.
    # CWL conformance uses alpine images that do not have bash.
    java \
        -Xmx2g \
        -Dconfig.file="${CROMWELL_BUILD_CROMWELL_CONFIG}" \
        -Dcall-caching.enabled=false \
        -Dsystem.job-shell=/bin/sh \
        -jar "${CROMWELL_BUILD_CROMWELL_JAR}" \
        server &

    CROMWELL_BUILD_CONFORMANCE_CROMWELL_PID=$!

    popd > /dev/null

    cromwell::private::add_exit_function cromwell::private::kill_conformance_cromwell
}

cromwell::private::kill_conformance_cromwell() {
    if [[ -n "${CROMWELL_BUILD_CONFORMANCE_CROMWELL_PID+set}" ]]; then
        cromwell::build::kill_tree "${CROMWELL_BUILD_CONFORMANCE_CROMWELL_PID}"
    fi
}

cromwell::private::run_conformance_wdl() {
    pushd "${CROMWELL_BUILD_CWL_TEST_RESOURCES}" > /dev/null

    java \
        -Xmx6g \
        -Dbackend.providers.Local.config.concurrent-job-limit="${CROMWELL_BUILD_CWL_TEST_PARALLELISM}" \
        -jar "${CROMWELL_BUILD_CROMWELL_JAR}" \
        run "${CROMWELL_BUILD_CWL_TEST_WDL}" \
        -i "${CROMWELL_BUILD_CWL_TEST_INPUTS}"

    popd > /dev/null
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
    cromwell::private::create_database_variables
    cromwell::private::verify_secure_build
    cromwell::private::make_build_directories
    cromwell::private::install_git_secrets
    cromwell::private::install_minnie_kenny
    cromwell::private::install_wait_for_it
    cromwell::private::setup_secure_resources

    case "${CROMWELL_BUILD_PROVIDER}" in
        "${CROMWELL_BUILD_PROVIDER_TRAVIS}")
            cromwell::private::stop_travis_defaults
            cromwell::private::delete_boto_config
            cromwell::private::delete_sbt_boot
            cromwell::private::upgrade_pip
            cromwell::private::pull_common_docker_images
            cromwell::private::start_docker_databases
            ;;
        "${CROMWELL_BUILD_PROVIDER_JENKINS}")
            cromwell::private::delete_boto_config
            cromwell::private::delete_sbt_boot
            cromwell::private::upgrade_pip
            ;;
        *)
            cromwell::private::pull_common_docker_images
            ;;
    esac
}

cromwell::build::setup_centaur_environment() {
    cromwell::private::create_centaur_variables
    cromwell::private::start_build_heartbeat
    cromwell::private::start_cromwell_log_tail
    cromwell::private::start_centaur_log_tail
    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        cromwell::private::add_exit_function cromwell::private::cat_centaur_log
    fi
}

cromwell::build::setup_conformance_environment() {
    cromwell::private::create_centaur_variables
    cromwell::private::create_conformance_variables
    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]]; then
        cromwell::private::install_cwltest
    fi
    cromwell::private::checkout_pinned_cwl
    cromwell::private::write_cwl_test_inputs
    cromwell::private::start_build_heartbeat
    cromwell::private::add_exit_function cromwell::private::cat_conformance_log
}

cromwell::build::setup_docker_environment() {
    cromwell::private::start_build_heartbeat

    if [[ "${CROMWELL_BUILD_PROVIDER}" == "${CROMWELL_BUILD_PROVIDER_TRAVIS}" ]]; then
        # Upgrade docker-compose so that we get the correct exit codes
        docker-compose -version
        sudo rm /usr/local/bin/docker-compose
        curl \
            -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" \
            > docker-compose
        chmod +x docker-compose
        sudo mv docker-compose /usr/local/bin
        docker-compose -version
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

cromwell::build::build_cromwell_docker() {
    if [[ "${CROMWELL_BUILD_IS_CI}" == "true" ]] || ! cromwell::private::exists_cromwell_docker; then
        echo "Please wait, building cromwell docker…"
        cromwell::private::build_cromwell_docker
    fi
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
    # Handle empty arrays in older versions of bash
    # https://stackoverflow.com/questions/7577052/bash-empty-array-expansion-with-set-u#answer-7577209
    "${CROMWELL_BUILD_ROOT_DIRECTORY}/centaur/test_cromwell.sh" \
        -n "${CROMWELL_BUILD_CENTAUR_CONFIG}" \
        -l "${CROMWELL_BUILD_LOG_DIRECTORY}" \
        -g \
        ${additional_args[@]+"${additional_args[@]}"} \
        "$@"
}

cromwell::build::run_conformance() {
    cromwell::private::start_conformance_cromwell

    # Give cromwell time to start up
    sleep 30

    cromwell::private::run_conformance_wdl
}

cromwell::build::generate_code_coverage() {
    if [[ "${CROMWELL_BUILD_GENERATE_COVERAGE}" == "true" ]]; then
        cromwell::private::generate_code_coverage
    fi
}

cromwell::build::publish_artifacts() {
    if [[ "${CROMWELL_BUILD_PROVIDER}" == "${CROMWELL_BUILD_PROVIDER_TRAVIS}" ]] && \
        [[ "${CROMWELL_BUILD_TYPE}" == "sbt" ]] && \
        [[ "${CROMWELL_BUILD_EVENT}" == "push" ]]; then

        if [[ "${CROMWELL_BUILD_BRANCH}" == "develop" ]]; then
            # Publish images for both the "cromwell develop branch" and the "cromwell dev environment".
            CROMWELL_SBT_DOCKER_TAGS=develop,dev \
                cromwell::private::publish_artifacts_and_docker \
                -Dproject.isSnapshot=true
            cromwell::private::push_publish_complete

        elif [[ "${CROMWELL_BUILD_BRANCH}" =~ ^[0-9\.]+_hotfix$ ]]; then
            # Docker tags float. "30" is the latest hotfix. Those dockers are published here on each hotfix commit.
            cromwell::private::publish_artifacts_and_docker -Dproject.isSnapshot=false

        elif [[ -n "${CROMWELL_BUILD_TAG:+set}" ]]; then
            # Artifact tags are static. Once "30" is set that is only "30" forever. Those artifacts are published here.
            cromwell::private::publish_artifacts_only \
                -Dproject.version="${CROMWELL_BUILD_TAG}" \
                -Dproject.isSnapshot=false

        elif [[ "${CROMWELL_BUILD_IS_SECURE}" == "true" ]]; then
            cromwell::private::publish_artifacts_check

        fi

    fi
}

cromwell::build::exec_retry_function() {
    local retried_function
    local retry_count
    local attempt
    local exit_status

    retried_function="${1:?exec_retry_function called without a function}"; shift
    retry_count="${1:-3}"; shift || true
    sleep_seconds="${1:-15}"; shift || true

    # https://unix.stackexchange.com/a/82610
    # https://stackoverflow.com/a/17336953
    for attempt in $(seq 0 "${retry_count}"); do
        [[ ${attempt} -gt 0 ]] && sleep "${sleep_seconds}"
        ${retried_function} && exit_status=0 && break || exit_status=$?
    done
    return ${exit_status}
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

cromwell::build::start_build_heartbeat() {
    cromwell::private::start_build_heartbeat
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
