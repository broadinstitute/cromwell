#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

mkdir -p mnt

# Read the service account credentials from vault:
read_service_account_from_vault

# Start up the CloudSQL DB. There is an --async on `patch` but it doesn't return a waitable ID...
gcloud_run_as_service_account "perf_sql_start_gcloud_${BUILD_NUMBER}" "gcloud -q --project broad-dsde-cromwell-perf sql instances patch ${CLOUD_SQL_INSTANCE} --activation-policy ALWAYS"
