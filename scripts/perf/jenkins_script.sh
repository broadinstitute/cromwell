#!/bin/bash
set -e

GCS_BUCKET=gs://cromwell-perf-test/

VAULT_TOKEN=$(cat /etc/vault-token-dsde)

DOCKER_ETC_PATH=/usr/share/etc

mkdir -p mnt

#TODO: need this in the google cloud docker in order to auth
#get startup script to run on new instance
curl https://raw.githubusercontent.com/broadinstitute/cromwell/$CROMWELL_PERF_SCRIPTS_BRANCH/scripts/perf/startup_script.sh -o mnt/startup_script.sh

DB_PASS=`docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN \
	broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf | jq '.data.db_pass'`

docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf/service-account-deployer | jq -r '.data.service_account' > mnt/sa.json
	
function join() { local IFS=","; echo "$*"; }

CLOUD_SQL_INSTANCE=cromwell-db-perf-test-$BUILD_TAG

# Clone the CloudSQL DB
# Note: Cloning the same database in parallel doesn't work.
# By doing it here we can run the jenkins jobs sequentially and ensure the database is cloned once at a time as well
docker run --name perf_sql_create_gcloud_$BUILD_NUMBER -v "$(pwd)"/mnt:$DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/sa.json &&\
    gcloud --project broad-dsde-cromwell-perf sql instances clone cromwell-perf-testing-base-09-24-18 ${CLOUD_SQL_INSTANCE}"
    
# Add a cromwell user to the cloned database
docker run --name perf_sql_user_gcloud_$BUILD_NUMBER -v "$(pwd)"/mnt:$DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/sa.json &&\
    gcloud --project broad-dsde-cromwell-perf sql users create cromwell --instance=${CLOUD_SQL_INSTANCE} --password=${DB_PASS}"

metadata=(
  "BUILD_NUMBER=$BUILD_NUMBER"
  "BUILD_TAG=$BUILD_TAG"
  "CLEAN_UP=$CLEAN_UP"
  "CLOUD_SQL_INSTANCE=$CLOUD_SQL_INSTANCE"
  "CROMWELL_DB_USER=cromwell"
  "CROMWELL_DB_PASS=$DB_PASS"
  "CROMWELL_DOCKER_IMAGE=$CROMWELL_DOCKER_IMAGE"
  "CROMWELL_PROJECT=broad-dsde-cromwell-perf"
  "CROMWELL_BUCKET=$GCS_BUCKET"
  "CROMWELL_STATSD_HOST=10.128.0.4"
  "CROMWELL_STATSD_PORT=8125"
  "CROMWELL_PERF_SCRIPTS_BRANCH=$CROMWELL_PERF_SCRIPTS_BRANCH"
  "GCS_REPORT_BUCKET=$GCS_REPORT_BUCKET"
  "GCS_REPORT_PATH=$GCS_REPORT_PATH"
  "TEST_CASE_DIRECTORY=$TEST_CASE_DIRECTORY"
)

docker run --name perf_gcloud_$BUILD_NUMBER -v "$(pwd)"/mnt:$DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/sa.json &&\
gcloud \
    --verbosity info \
    --project broad-dsde-cromwell-perf \
    compute \
    instances \
    create $BUILD_TAG \
    --zone us-central1-c \
    --source-instance-template cromwell-perf-template-09-18 \
    --metadata-from-file startup-script=$DOCKER_ETC_PATH/startup_script.sh \
    --metadata \
        $(join ${metadata[@]})"
