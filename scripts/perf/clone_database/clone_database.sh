#!/bin/bash
set -e

VAULT_TOKEN=$(cat /etc/vault-token-dsde)

DOCKER_ETC_PATH=/usr/share/etc

mkdir -p mnt

# Read the DB password from vault:
DB_PASS=`docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN \
	broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf | jq '.data.db_pass'`

# Read the service account credentials from vault:
docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf/service-account-deployer | jq -r '.data.service_account' > mnt/sa.json

# Clone the CloudSQL DB
# Note: Cloning the same database in parallel doesn't work.
# By doing it here we can run the jenkins jobs sequentially and ensure the database is cloned once at a time as well
SQL_OPERATION=$(docker run --name perf_sql_create_gcloud_$BUILD_NUMBER -v "$(pwd)"/mnt:$DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/sa.json 2> /dev/null &&\
    gcloud --project broad-dsde-cromwell-perf sql instances clone --async ${CLOUD_SQL_INSTANCE_TO_CLONE} ${CLOUD_SQL_INSTANCE_NEW_NAME} --format='value(name)'")

# The cloning itself sometimes takes a long time and the clone command errors out when that happens.
# Instead use the --async flag in the clone command above and then explicitly wait for the operation to be done. Timeout 15 minutes
docker run --name perf_sql_create_gcloud_$BUILD_NUMBER -v "$(pwd)"/mnt:$DOCKER_ETC_PATH -e DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/sa.json &&\
    gcloud beta sql operations wait --timeout=900 --project broad-dsde-cromwell-perf $SQL_OPERATION"

# Add a cromwell user to the cloned database
docker run --name perf_sql_user_gcloud_$BUILD_NUMBER -v "$(pwd)"/mnt:$DOCKER_ETC_PATH --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/sa.json &&\
    gcloud --project broad-dsde-cromwell-perf sql users create cromwell --instance=${CLOUD_SQL_INSTANCE_NEW_NAME} --password=${DB_PASS}"

echo ""
