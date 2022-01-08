#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

mkdir -p mnt

# Read the service account credentials from vault:
read_service_account_from_vault

gcloud_run_as_service_account "perf_sql_delete_gcloud_${BUILD_NUMBER}" " \
  gcloud \
      --verbosity info \
      --project broad-dsde-cromwell-perf \
      sql instances delete ${CLOUD_SQL_INSTANCE} -q"
