#!/usr/bin/env bash

set -x

wait_for_cromwell() {
  git clone https://github.com/vishnubob/wait-for-it.git /wait-for-it
  chmod u+x /wait-for-it/wait-for-it.sh
  # Give 1 minute to cromwell to be online
  /wait-for-it/wait-for-it.sh -h localhost -p 8000 -t 60

  # TODO install on the image
  apt-get update && apt-get install --assume-yes jq
  
  CROMWELL_VERSION=$(curl -X GET "http://localhost:8000/engine/v1/version" -H "accept: application/json" | jq -r '.cromwell')
  if [ -z ${CROMWELL_VERSION} ]
  then
    CROMWELL_VERSION="unknown"
  fi
  export CROMWELL_VERSION
}

export_logs() {
    # Copy the cromwell container logs
    gsutil -h "Content-Type:text/plain" cp $(docker inspect --format='{{.LogPath}}' cromwell) "${REPORT_PATH}/cromwell.log"
    # Copy the docker daemon logs
    gsutil -h "Content-Type:text/plain" cp /var/log/daemon.log "${REPORT_PATH}/daemon.log"
}

clean_up() {
    gcloud sql instances delete cromwell-db-${BUILD_ID}
    gcloud compute instances delete $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=us-central1-c -q
}

run_test() {
    export CENTAUR_TEST_NAME=$(grep name: ${CENTAUR_TEST_FILE} | cut -d ' ' -f 2)
    export REPORT_PATH="${CENTAUR_TEST_NAME}/${CROMWELL_VERSION}/${BUILD_ID}"
    export REPORT_URL="gs://${REPORT_BUCKET}/${REPORT_PATH}"

    # Ideally this image is in sync with CROMWELL_BRANCH so that there's nothing to (re-)compile here
    gcloud docker -- pull us.gcr.io/broad-dsde-cromwell-perf/centaur:perf
    
    docker run \
     -v ${CROMWELL_ROOT}:${CROMWELL_ROOT} \
     -e CENTAUR_TEST_FILE \
     -e REPORT_BUCKET \
     -e REPORT_PATH \
     -e CROMWELL_BRANCH \
     -e PERF_ROOT \
     --network=vm_scripts_mysql_net \
     us.gcr.io/broad-dsde-cromwell-perf/centaur:perf /bin/bash /app/scripts/perf/vm_scripts/centaur/test_runner.sh
}
