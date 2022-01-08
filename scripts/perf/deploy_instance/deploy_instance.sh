#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

DB_PASS=$(read_path_from_vault_json "secret/dsp/cromwell/perf" '.data.db_pass')

read_service_account_from_vault
	
function join() { local IFS=","; echo "$*"; }

metadata=(
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
)

cp scripts/perf/deploy_instance/run_on_instance.sh mnt/

gcloud_run_as_service_account "perf_deploy_instance_${BUILD_NUMBER}" \
  "gcloud \
    --verbosity info \
    --project broad-dsde-cromwell-perf \
    compute \
    instances \
    create ${CROMWELL_INSTANCE_NAME} \
    --zone us-central1-c \
    --source-instance-template $INSTANCE_TEMPLATE \
    --metadata-from-file startup-script=$DOCKER_ETC_PATH/run_on_instance.sh \
    --metadata \
        $(join ${metadata[@]})" | tee dockerResult.txt

typeset CROMWELL_UNDER_TEST=$(cat dockerResult.txt | tail -n1 | awk '{print $5}' )
typeset CROMWELL_VALUE_TO_RETURN=$(cat dockerResult.txt | tail -n1 | awk "{print \$${FIELD_NO_TO_RETURN}}" )

if test -z "CROMWELL_UNDER_TEST" -o -z "CROMWELL_VALUE_TO_RETURN"
then
  echo "One of CROMWELL_UNDER_TEST or CROMWELL_VALUE_TO_RETURN are empty ('${CROMWELL_UNDER_TEST}', '${CROMWELL_VALUE_TO_RETURN}')"
  exit 1
else
  echo "Determined that CROMWELL_UNDER_TEST=${CROMWELL_UNDER_TEST}, CROMWELL_VALUE_TO_RETURN=${CROMWELL_VALUE_TO_RETURN}"
fi

custom_wait_for_cromwell "${CROMWELL_UNDER_TEST}"

mkdir -p output
echo "CROMWELL_UNDER_TEST=${CROMWELL_VALUE_TO_RETURN}" > output/cromwell.properties
