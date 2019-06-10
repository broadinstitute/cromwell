#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

GCS_BUCKET=gs://cromwell-perf-test/

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
  "CROMWELL_BUCKET=$GCS_BUCKET"
  "CROMWELL_STATSD_HOST=10.128.0.4"
  "CROMWELL_STATSD_PORT=8125"
  "CROMWELL_PERF_SCRIPTS_BRANCH=${REPO_BRANCH}"
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

if test -z "$CROMWELL_UNDER_TEST"
then
  echo "\$CROMWELL_UNDER_TEST is empty"
  exit 1
else
  echo "Determined that CROMWELL_UNDER_TEST=${CROMWELL_UNDER_TEST}"
fi

custom_wait_for_cromwell

mkdir -p output
echo "CROMWELL_UNDER_TEST=${CROMWELL_UNDER_TEST}" > output/cromwell.properties
