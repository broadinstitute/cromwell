#!/bin/bash
set -e
VAULT_TOKEN=$(cat /etc/vault-token-dsde)
DB_PASS=`docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN \
	broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf | jq '.data.db_pass'`

gcloud \
    compute \
    instances \
    create perf-test-$BUILD_NUMBER \
    --project broad-dsde-cromwell-perf \
    --zone us-central1-c \
    --source-instance-template cromwell-perf-template-update \
    --metadata-from-file startup-script=startup_script.sh \
    --metadata \
        CROMWELL_DB_USER=cromwell,CROMWELL_DB_PASS=$DB_PASS,CLOUD_SQL_INSTANCE=cromwell-perf-testing-db-$BUILD_NUMBER,CROMWELL_VERSION=34,CROMWELL_PROJECT=broad-dsde-cromwell-perf,CROMWELL_BUCKET=gs://debtest3/,CROMWELL_STATSD_HOST=broad.io/batch-grafana,CROMWELL_STATSD_PORT=8125
