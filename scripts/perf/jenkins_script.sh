#!/bin/bash
set -e
VAULT_TOKEN=$(cat /etc/vault-token-dsde)
DB_PASS=docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN \
	broadinstitute/dsde-toolbox vault read -format=json secret/dsp/cromwell/perf | jq '.db_pass'

gcloud --project broad-dsde-cromwell-perf sql instances clone cromwell-perf-testing-base cromwell-perf-testing-$BUILD_NUMBER
gcloud --project broad-dsde-cromwell-perf sql users create cromwell --instance=cromwell-perf-testing-db-$BUILD_NUMBER --password=$DB_PASS
gcloud \
    --project broad-dsde-cromwell-perf \
    compute \
    instances \
    create perf-test-$BUILD_NUMBER \
    --source-instance-template cromwell-perf-template-update \
    --metadata-from-file startup-script=startup_script.sh \
    --metadata \
        CROMWELL_DB_USER=cromwell,CROMWELL_DB_PASS=$DB_PASS,CLOUD_SQL_INSTANCE=cromwell-perf-testing-db-$BUILD_NUMBER,CROMWELL_VERSION=34,CROMWELL_PROJECT=broad-dsde-cromwell-perf,CROMWELL_BUCKET=gs://debtest3/,CROMWELL_STATSD_HOST=broad.io/batch-grafana,CROMWELL_STATSD_PORT=8125
