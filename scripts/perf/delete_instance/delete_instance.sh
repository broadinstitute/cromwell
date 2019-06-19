#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

mkdir -p mnt

# Read the service account credentials from vault:
read_service_account_from_vault

gcloud_run_as_service_account "perf_delete_gcloud_instance_${BUILD_NUMBER}" " \
  gcloud \
    --verbosity info \
    --project broad-dsde-cromwell-perf \
    compute instances delete ${CROMWELL_INSTANCE_NAME} --zone=us-central1-c -q"
