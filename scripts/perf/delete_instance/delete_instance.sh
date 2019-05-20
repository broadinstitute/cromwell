#!/bin/bash
set -e

VAULT_TOKEN=$(cat /etc/vault-token-dsde)

DOCKER_ETC_PATH=/usr/share/etc

mkdir -p mnt

# Read the service account credentials from vault:
docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf/service-account-deployer | jq -r '.data.service_account' > mnt/sa.json

docker run --name perf_sql_delete_gcloud_${BUILD_NUMBER} \
  -v "$(pwd)"/mnt:$DOCKER_ETC_PATH \
  -e DOCKER_ETC_PATH="${DOCKER_ETC_PATH}" \
  -e CROMWELL_INSTANCE_NAME="${CROMWELL_INSTANCE_NAME}" \
  --rm \
  google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file ${DOCKER_ETC_PATH}/sa.json && \
    gcloud \
      --verbosity info \
      --project broad-dsde-cromwell-perf \
      compute instances delete ${CROMWELL_INSTANCE_NAME} --zone=us-central1-c -q"
