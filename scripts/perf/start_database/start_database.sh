#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

mkdir -p mnt

# Read the service account credentials from vault:
read_service_account_from_vault

# Start up the CloudSQL DB
SQL_OPERATION=$(gcloud_run_as_service_account "perf_sql_start_gcloud_${BUILD_NUMBER}" "gcloud -q --project broad-dsde-cromwell-perf sql instances patch --async ${CLOUD_SQL_INSTANCE_TO_START} --activation-policy ALWAYS")

gcloud_run_as_service_account "perf_sql_await_gcloud_${BUILD_NUMBER}" "gcloud beta sql operations wait --timeout=900 --project broad-dsde-cromwell-perf ${SQL_OPERATION}"
