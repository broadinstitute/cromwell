#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

DB_PASS=$(read_path_from_vault_json "secret/dsp/cromwell/perf" '.data.db_pass')

read_service_account_from_vault

runner_metadata=(
  "BUILD_NUMBER=$BUILD_NUMBER"
  "BUILD_TAG=$BUILD_TAG"
  "CLEAN_UP=false"
  "CLOUD_SQL_INSTANCE=$CLOUD_SQL_INSTANCE"
  "CROMWELL_DB_USER=cromwell"
  "CROMWELL_DB_PASS=$DB_PASS"
  "CROMWELL_DOCKER_IMAGE=$CROMWELL_DOCKER_IMAGE"
  "CROMWELL_PROJECT=broad-dsde-cromwell-perf"
  "CROMWELL_BUCKET=${CROMWELL_EXECUTION_BUCKET}"
  "CROMWELL_STATSD_HOST=10.128.0.4"
  "CROMWELL_STATSD_PORT=8125"
  "CROMWELL_PERF_SCRIPTS_BRANCH=${REPO_BRANCH}"
  "GCS_REPORT_PATH=${GCS_REPORT_PATH}"
  "GCS_REPORT_BUCKET=${GCS_REPORT_BUCKET}"
  "TEST_CASE_DIRECTORY=${TEST_CASE_DIRECTORY}"
  "CONFIG_OVERRIDE"="-Dservices.MetadataService.config.metadata-summary-refresh-interval = \"Inf\""
)

reader_metadata=(
  "BUILD_NUMBER=$BUILD_NUMBER"
  "BUILD_TAG=$BUILD_TAG"
  "CLEAN_UP=false"
  "CLOUD_SQL_INSTANCE=$CLOUD_SQL_INSTANCE"
  "CROMWELL_DB_USER=cromwell"
  "CROMWELL_DB_PASS=$DB_PASS"
  "CROMWELL_DOCKER_IMAGE=$CROMWELL_DOCKER_IMAGE"
  "CROMWELL_PROJECT=broad-dsde-cromwell-perf"
  "CROMWELL_BUCKET=${CROMWELL_EXECUTION_BUCKET}"
  "CROMWELL_STATSD_HOST=10.128.0.4"
  "CROMWELL_STATSD_PORT=8125"
  "CROMWELL_PERF_SCRIPTS_BRANCH=${REPO_BRANCH}"
  "GCS_REPORT_PATH=${GCS_REPORT_PATH}"
  "GCS_REPORT_BUCKET=${GCS_REPORT_BUCKET}"
  "TEST_CASE_DIRECTORY=reader"
  "CONFIG_OVERRIDE"=""
)

cp scripts/perf/deploy_instance/run_on_instance.sh mnt/

gcloud_deploy_instance \
  "perf_deploy_instance_${BUILD_NUMBER}-runner101" \
   "${CROMWELL_INSTANCE_NAME}-runner101" \
    "${INSTANCE_TEMPLATE}" \
     "$(join ${runner_metadata[@]})" | tee runner1.txt

gcloud_deploy_instance \
  "perf_deploy_instance_${BUILD_NUMBER}-runner102" \
   "${CROMWELL_INSTANCE_NAME}-runner102" \
   "${INSTANCE_TEMPLATE}" \
   "$(join ${runner_metadata[@]})" | tee runner2.txt

# Note: The reader_metadata currently allows the default summarizer to continue running...
# ... that's fine as long as there's only one reader, but if we add more, we should make the "summarizer" a separate VM.

gcloud_deploy_instance \
  "perf_deploy_instance_${BUILD_NUMBER}-reader101" \
   "${CROMWELL_INSTANCE_NAME}-reader101" \
   "${INSTANCE_TEMPLATE}" \
   "$(join ${reader_metadata[@]})" | tee reader1.txt

typeset RUNNER_UNDER_TEST_1=$(cat runner1.txt | tail -n1 | awk '{print $5}')
typeset RUNNER_UNDER_TEST_2=$(cat runner2.txt | tail -n1 | awk '{print $5}')

typeset READER_UNDER_TEST_1=$(cat reader1.txt | tail -n1 | awk '{print $5}')
typeset CROMWELL_VALUE_TO_RETURN=$(cat reader1.txt | tail -n1 | awk "{print \$${FIELD_NO_TO_RETURN}}" )

if test -z "${RUNNER_UNDER_TEST_1}" -o -z "${RUNNER_UNDER_TEST_2}" -o -z "${READER_UNDER_TEST_1}"
then
  echo "One of RUNNER_UNDER_TEST_1 or RUNNER_UNDER_TEST_2 or READER_UNDER_TEST_1 are empty ('${RUNNER_UNDER_TEST_1}', '${RUNNER_UNDER_TEST_2}, ${READER_UNDER_TEST_1}')"
  exit 1
else
  echo "Determined that RUNNER_UNDER_TEST_1=${RUNNER_UNDER_TEST_1}, RUNNER_UNDER_TEST_2=${RUNNER_UNDER_TEST_2}, READER_UNDER_TEST_1=${READER_UNDER_TEST_1}"
fi

custom_wait_for_cromwell "${RUNNER_UNDER_TEST_1}"
custom_wait_for_cromwell "${RUNNER_UNDER_TEST_2}"
custom_wait_for_cromwell "${READER_UNDER_TEST_1}"

mkdir -p output
echo "CROMWELL_UNDER_TEST=${CROMWELL_VALUE_TO_RETURN}" > output/cromwell.properties
