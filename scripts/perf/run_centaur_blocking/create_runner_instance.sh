#!/bin/bash

set -euo pipefail

source scripts/perf/helper.inc.sh

# Prepare the startup script to run on new instance
mkdir -p mnt
cp scripts/perf/run_centaur_blocking/test_runner_script.sh mnt/


read_service_account_from_vault

metadata=(
  "BUILD_NUMBER=${BUILD_NUMBER}"
  "BUILD_TAG=${BUILD_TAG}"
  "CLEAN_UP=true"
  "CROMWELL_PROJECT=${CROMWELL_PERF_PROJECT}"
  "CROMWELL_BUCKET=${CROMWELL_EXECUTION_BUCKET}"
  "CROMWELL_PERF_SCRIPTS_BRANCH=${CROMWELL_BRANCH}"
  "CROMWELL_UNDER_TEST=${CROMWELL_UNDER_TEST}"
  "GCS_REPORT_PATH=${GCS_REPORT_PATH}"
  "GCS_REPORT_BUCKET=${GCS_REPORT_BUCKET}"
  "TEST_CASE_DIRECTORY=${TEST_CASE_DIRECTORY}"
)

# NB: The 'cromwell-perf' tag is required to give cromwell the 'perf' firewall rules

gcloud_run_as_service_account "perf_centaur_test_runner_${BUILD_NUMBER}" \
  "gcloud \
    --verbosity info \
    --project ${CROMWELL_PERF_PROJECT} \
    compute \
    instances \
    create ${BUILD_TAG} \
    --zone us-central1-c \
    --source-instance-template ${INSTANCE_TEMPLATE} \
    --tags cromwell-perf \
    --metadata-from-file startup-script=${DOCKER_ETC_PATH}/test_runner_script.sh \
    --metadata \
        $(join ${metadata[@]})"  | tee startupResult.txt

typeset CENTAUR_TEST_RUNNER=$(cat startupResult.txt | tail -n1 | awk '{print $5}' )

echo "Determined that CENTAUR_TEST_RUNNER=${CENTAUR_TEST_RUNNER}"

wait_for_ping_to_start_then_stop ${CENTAUR_TEST_RUNNER}
