#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

VAULT_TOKEN=$(cat /etc/vault-token-dsde)

typeset DOCKER_ETC_PATH=/usr/share/etc

mkdir -p mnt

# Read the DB password from vault:
DB_PASS=$( docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN \
	broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf \
	| jq --exit-status '.data.db_pass' )

# Read the service account credentials from vault:
docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf/service-account-deployer | jq --exit-status -r '.data.service_account' > mnt/sa.json

# Clone the CloudSQL DB
# Note: Cloning the same database in parallel doesn't work.
# By doing it here we can run the jenkins jobs sequentially and ensure the database is cloned once at a time as well
SQL_OPERATION=$(gcloud_run_as_service_account "perf_sql_create_gcloud_${BUILD_NUMBER}" "gcloud --project broad-dsde-cromwell-perf sql instances clone --async ${CLOUD_SQL_INSTANCE_TO_CLONE} ${CLOUD_SQL_INSTANCE_NEW_NAME} --format='value(name)'" )

gcloud_run_as_service_account "perf_sql_await_gcloud_${BUILD_NUMBER}" "gcloud beta sql operations wait --timeout=900 --project broad-dsde-cromwell-perf ${SQL_OPERATION}"

gcloud_run_as_service_account "perf_sql_user_gcloud_${BUILD_NUMBER}" "gcloud --project broad-dsde-cromwell-perf sql users create cromwell --instance=${CLOUD_SQL_INSTANCE_NEW_NAME} --password=${DB_PASS}"
