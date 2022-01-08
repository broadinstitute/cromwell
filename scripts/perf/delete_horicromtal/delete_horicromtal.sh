#!/bin/bash
set -euo pipefail

source scripts/perf/helper.inc.sh

mkdir -p mnt

# Read the service account credentials from vault:
read_service_account_from_vault

for instance in runner101 runner102 reader101
do
    gcloud_delete_instance \
        "perf_delete_gcloud_instance_${BUILD_NUMBER}-${instance}" \
        "${CROMWELL_INSTANCE_NAME}-${instance}"
done
