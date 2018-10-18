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
    # Copy the cromwell container logs
    gsutil -h "Content-Type:text/plain" cp $(docker inspect --format='{{.LogPath}}' cromwell) "${REPORT_URL}/cromwell.log" || true
    # Copy the docker daemon logs
    gsutil -h "Content-Type:text/plain" cp /var/log/daemon.log "${REPORT_URL}/daemon.log" || true
    # Copy the statsd logs
    gsutil -h "Content-Type:text/plain" cp "${STATSD_PROXY_APP_DIR}/statsd.log" "${REPORT_URL}/statsd.log" || true
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

set_up() {
    apt-get update
    apt-get install --assume-yes apache2-utils
}

clean_up() {
    gcloud sql instances delete cromwell-db-${BUILD_ID}
    gcloud compute instances delete $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=us-central1-c -q
}

run_test() {
    cd ${CROMWELL_ROOT}

    export CENTAUR_TEST_NAME=$(grep name: ${CENTAUR_TEST_FILE} | cut -d ' ' -f 2)
    export REPORT_PATH="${CENTAUR_TEST_NAME}/${CROMWELL_VERSION}/${BUILD_ID}"
    export REPORT_URL="gs://${REPORT_BUCKET}/${REPORT_PATH}"
    
    sbt -Dconfig.file=${CENTAUR_CONF_DIR}/centaur.conf "centaur/it:testOnly centaur.ExternalTestCaseSpec"
}

shutdown() {
    export_logs
    docker-compose -f ${PERF_ROOT}/vm_scripts/docker-compose.yml down
    clean_up
}
