gcloud \
    --project broad-dsde-cromwell-perf \
    compute \
    instances \
    create perf-test \
    --source-instance-template cromwell-perf-template-update \
    --metadata-from-file startup-script=startup_script.sh \
    --metadata \
        CROMWELL_DB_USER=cromwell,CROMWELL_DB_PASSWORD=cromwell,CLOUD_SQL_INSTANCE=cromwell-perf-testing-db,CROMWELL_VERSION=34,CROMWELL_PROJECT=broad-dsde-cromwell-perf,CROMWELL_BUCKET=gs://debtest3/,CROMWELL_STATSD_HOST=broad.io/batch-grafana,CROMWELL_STATSD_PORT=8125
