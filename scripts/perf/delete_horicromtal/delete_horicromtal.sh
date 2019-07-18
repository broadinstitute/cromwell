#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

mkdir -p mnt

# Read the service account credentials from vault:
read_service_account_from_vault

gcloud_delete_instance \
  "perf_delete_gcloud_instance_${BUILD_NUMBER}-runner101" \
  "${CROMWELL_INSTANCE_NAME}-runner101"

gcloud_delete_instance \
  "perf_delete_gcloud_instance_${BUILD_NUMBER}-runner102" \
  "${CROMWELL_INSTANCE_NAME}-runner102"

gcloud_delete_instance \
  "perf_delete_gcloud_instance_${BUILD_NUMBER}-reader101" \
  "${CROMWELL_INSTANCE_NAME}-reader101"
