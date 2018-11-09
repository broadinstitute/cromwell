#!/usr/bin/env bash

set -x

wait_for_cromwell() {
  git clone https://github.com/vishnubob/wait-for-it.git /wait-for-it
  chmod u+x /wait-for-it/wait-for-it.sh
  # Give 5 minutes to cromwell to be online
  /wait-for-it/wait-for-it.sh localhost:8000 -t 300
  READY=$?
  if [ ${READY} -eq 0 ]
  then
      # Just wait a bit longer - no rush, this is a chill VM - this is because cromwell responds to requests before being really ready...
      sleep 30
    
      CROMWELL_VERSION=$(curl -X GET "http://localhost:8000/engine/v1/version" -H "accept: application/json" | jq -r '.cromwell')
      if [ -z ${CROMWELL_VERSION} ]
      then
        echo "Cromwell was up but failed to return its version, so something went wrong, shutting down"
        shutdown
      fi
      export CROMWELL_VERSION
  else
    echo "Cromwell was not ready after 5 minutes, shutting down"
    shutdown
  fi
}

export_logs() {
    export REPORT_URL="gs://${GCS_REPORT_BUCKET}/${GCS_REPORT_PATH}"

    # Copy the cromwell container logs
    gsutil -h "Content-Type:text/plain" cp $(docker inspect --format='{{.LogPath}}' cromwell) "${REPORT_URL}/cromwell.log" || true
    # Copy the docker daemon logs
    gsutil -h "Content-Type:text/plain" cp /var/log/daemon.log "${REPORT_URL}/daemon.log" || true
    # Copy the statsd logs
    gsutil -h "Content-Type:text/plain" cp "${STATSD_PROXY_APP_DIR}/statsd.log" "${REPORT_URL}/statsd.log" || true
    # Copy centaur log
    gsutil -h "Content-Type:text/plain" cp "${CROMWELL_ROOT}/centaur.log" "${REPORT_URL}/centaur.log" || true
    # Copy test_rc.log
    echo $TEST_RC > test_rc.txt && gsutil -h "Content-Type:text/plain" cp test_rc.txt "${REPORT_URL}/test_rc.txt" || true
}

prepare_statsd_proxy() {
    # Build the image
    cd ${CROMWELL_ROOT}
    export CROMWELL_SBT_DOCKER_TAGS=perf
    sbt statsDProxy/docker
    
    # Export configuration variables
    export STATSD_PROXY_APP_DIR="${PERF_ROOT}/vm_scripts/statsd-proxy"
    # Hard code those for now as there's no obvious reason to make them configurable
    export PROXY_HOST=0.0.0.0
    export PROXY_PORT=9125
}

clean_up() {
    if [ "${CLEAN_UP}" = true ]
    then
        gcloud sql instances delete ${CLOUD_SQL_INSTANCE}
        gcloud compute instances delete $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=us-central1-c -q
    fi
}

run_test() {
    cd ${CROMWELL_ROOT}

    if [ -z ${GCS_REPORT_PATH} ]
    then
        export GCS_REPORT_PATH="${TEST_CASE_DIRECTORY}/${CROMWELL_VERSION}/${BUILD_ID}"
    fi

    sbt -Dconfig.file=${CENTAUR_CONF_DIR}/centaur.conf "centaur/it:testOnly centaur.ExternalTestCaseSpec" | tee centaur.log
    export TEST_RC=$?
}

shutdown() {
    export_logs
    docker-compose -f ${PERF_ROOT}/vm_scripts/docker-compose.yml down
    clean_up
}
