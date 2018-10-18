#!/bin/bash
set -e

GCS_BUCKET=gs://cromwell-perf-test/

VAULT_TOKEN=$(cat /etc/vault-token-dsde)

DOCKER_ETC_PATH=/usr/share/etc

mkdir -p mnt

#TODO: need this in the google cloud docker in order to auth
#get startup script to run on new instance
curl https://raw.githubusercontent.com/broadinstitute/cromwell/$CROMWELL_BRANCH/scripts/perf/startup_script.sh > mnt/startup_script.sh

DB_PASS=`docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN \
	broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf | jq '.data.db_pass'`

docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf/service-account-deployer | jq -r '.data.service_account' > mnt/sa.json

DB_PASS=`docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN \
	broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf | jq '.data.db_pass'`

docker run --name perf_gcloud_$BUILD_NUMBER -v "$(pwd)"/mnt:$DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/sa.json &&\
gcloud \
    --verbosity info \
    --project broad-dsde-cromwell-perf \
    compute \
    instances \
    create perf-test-$BUILD_NUMBER \
    --zone us-central1-c \
    --source-instance-template cromwell-perf-template-09-18 \
    --metadata-from-file startup-script=$DOCKER_ETC_PATH/startup_script.sh \
    --metadata \
        CROMWELL_DB_USER=cromwell,CROMWELL_DB_PASS=$DB_PASS,CLOUD_SQL_INSTANCE=cromwell-db-perf-test-$BUILD_NUMBER,DO_SHUTDOWN=$DO_SHUTDOWN,CROMWELL_IMAGE=$CROMWELL_VERSION_NUMBER,CROMWELL_PROJECT=broad-dsde-cromwell-perf,CROMWELL_BUCKET=$GCS_BUCKET,CROMWELL_STATSD_HOST=10.128.0.4,CROMWELL_STATSD_PORT=8125,CROMWELL_STATSD_PREFIX=perf-test-$BUILD_NUMBER,CROMWELL_BRANCH_NAME=$CROMWELL_BRANCH,WORKFLOW_SCRIPT=$WORKFLOW_SCRIPT"
