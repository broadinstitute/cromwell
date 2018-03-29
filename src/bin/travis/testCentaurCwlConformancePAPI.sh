#!/usr/bin/env bash

if [ "$TRAVIS_SECURE_ENV_VARS" = "false" ]; then
    echo "************************************************************************************************"
    echo "************************************************************************************************"
    echo "**                                                                                            **"
    echo "**  WARNING: Encrypted keys are unavailable to automatically test JES with centaur. Exiting.  **"
    echo "**                                                                                            **"
    echo "************************************************************************************************"
    echo "************************************************************************************************"
    exit 0
fi

set -e

sudo -H pip install --upgrade pip
sudo -H pip install cwltest

# -- BEGIN PAPI related conf generation

# TURN OFF LOGGING WHILE WE TALK TO DOCKER/VAULT
set +x

# Login to docker to access the dsde-toolbox
docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"

# Login to vault to access secrets
docker run --rm \
    -v "$HOME:/root:rw" \
    broadinstitute/dsde-toolbox \
    vault auth "$JES_TOKEN" < /dev/null > /dev/null && echo vault auth success

# For now, leave logging off
#set -x

# Render secrets
docker run --rm \
    -v "$HOME:/root:rw" \
    -v "$PWD/src/bin/travis/resources:/working" \
    -v "$PWD:/output" \
    -e ENVIRONMENT=not_used \
    -e INPUT_PATH=/working \
    -e OUT_PATH=/output \
    broadinstitute/dsde-toolbox render-templates.sh

JES_CONF="$(pwd)/jes_centaur.conf"

# -- END PAPI related conf generation

# Turn off call caching as hashing doesn't work since it sees local and not GCS paths.
sed -i '/^call-caching\s*/{N;s/enabled.*/  enabled: false/;}' ${JES_CONF}

ENABLE_COVERAGE=true sbt assembly

git clone https://github.com/common-workflow-language/common-workflow-language.git
cd common-workflow-language
# checkout a known git hash to prevent the tests from changing out from under us
git checkout 1646a398eacee14899d56945204c0162372fb5d8
cd ..

CROMWELL_JAR=$(find "$(pwd)/server/target/scala-2.12" -name "cromwell-*.jar")
CENTAUR_CWL_RUNNER="$(pwd)/centaurCwlRunner/src/bin/centaur-cwl-runner.bash"
CENTAUR_CWL_RUNNER_MODE="papi"
GOOGLE_AUTH_MODE="service-account"
GOOGLE_SERVICE_ACCOUNT_JSON="$(pwd)/cromwell-service-account.json"
CONFORMANCE_EXPECTED_FAILURES=$(pwd)/src/bin/travis/resources/papi_conformance_expected_failures.txt
CWL_DIR=$(pwd)/common-workflow-language
CWL_TEST_DIR=$(pwd)/common-workflow-language/v1.0/v1.0
CWL_CONFORMANCE_TEST_WDL=$(pwd)/src/bin/travis/resources/cwl_conformance_test.wdl
CWL_CONFORMANCE_TEST_INPUTS=$(pwd)/cwl_conformance_test.inputs.json
CWL_CONFORMANCE_TEST_PARALLELISM=10 # Set too high and it will cause false negatives due to cromwell server timeouts.
PAPI_INPUT_GCS_PREFIX=gs://centaur-cwl-conformance/cwl-inputs/

CURRENT_DIR=$(pwd)
# The java programs must run from the directory that contains a number of the cwl imports or the tests won't pass.
cd ${CWL_TEST_DIR}

# Export variables used in conf files and commands
export CENTAUR_CWL_RUNNER_MODE
export GOOGLE_AUTH_MODE
export GOOGLE_SERVICE_ACCOUNT_JSON
export PAPI_INPUT_GCS_PREFIX

java \
  -Xmx2g \
  -Dconfig.file="$JES_CONF" \
  -Dsystem.new-workflow-poll-rate=1 \
  -jar "${CROMWELL_JAR}" server &

CROMWELL_PID=$$

sleep 30

cat <<JSON >${CWL_CONFORMANCE_TEST_INPUTS}
{
  "cwl_conformance_test.cwl_dir": "$CWL_DIR",
  "cwl_conformance_test.centaur_cwl_runner": "$CENTAUR_CWL_RUNNER",
  "cwl_conformance_test.conformance_expected_failures": "$CONFORMANCE_EXPECTED_FAILURES"
}
JSON

java \
  -Xmx2g \
  -Dbackend.providers.Local.config.concurrent-job-limit=${CWL_CONFORMANCE_TEST_PARALLELISM} \
  -jar ${CROMWELL_JAR} \
  run ${CWL_CONFORMANCE_TEST_WDL} \
  -i ${CWL_CONFORMANCE_TEST_INPUTS}

cd $CURRENT_DIR

if [ "$TRAVIS_EVENT_TYPE" != "cron" ]; then
    sbt coverageReport --warn
    sbt coverageAggregate --warn
    bash <(curl -s https://codecov.io/bash) >/dev/null
fi
