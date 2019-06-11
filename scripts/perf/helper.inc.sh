#!/usr/bin/env bash

# Set some global constants:
typeset DOCKER_ETC_PATH="/usr/share/etc"

if test -f "/etc/vault-token-dsde"
then
  typeset VAULT_TOKEN=$(cat /etc/vault-token-dsde)
fi


wait_for_cromwell() {
  git clone https://github.com/vishnubob/wait-for-it.git /wait-for-it
  chmod u+x /wait-for-it/wait-for-it.sh
  # Give 5 minutes to cromwell to be online

  echo "[$(date)] Waiting for http://${CROMWELL_UNDER_TEST}:8000/engine/v1/version to appear..."
  /wait-for-it/wait-for-it.sh ${CROMWELL_UNDER_TEST}:8000 -t 300
  READY=$?
  if [ ${READY} -eq 0 ]
  then
      # Just wait a bit longer - no rush, this is a chill VM - this is because cromwell responds to requests before being really ready...
      sleep 30
    
      CROMWELL_VERSION=$(curl -X GET "http://${CROMWELL_UNDER_TEST}:8000/engine/v1/version" -H "accept: application/json" | jq -r '.cromwell')
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

# Waits for Cromwell. Works even from *outside* a google instance
custom_wait_for_cromwell() {

  set +e

  RESULT=1
  ATTEMPTS=0
  MAX_ATTEMPTS=20

  while [ "${ATTEMPTS}" -le "${MAX_ATTEMPTS}" -a "${RESULT}" -gt "0" ]
  do
    echo "[$(date)] Waiting for Cromwell (http://${CROMWELL_UNDER_TEST}:8000/engine/v1/version) to come up (tried ${ATTEMPTS} times so far)"
    sleep 30
    ATTEMPTS=$((ATTEMPTS + 1))

    CROMWELL_VERSION_JSON=$(curl -X GET "http://${CROMWELL_UNDER_TEST}:8000/engine/v1/version" -H "accept: application/json")
    RESULT=$?

    CROMWELL_VERSION=$(echo "${CROMWELL_VERSION_JSON}" | jq -r '.cromwell')

    if [ "${RESULT}" -eq "0" -a -z "${CROMWELL_VERSION}" ]
    then
      echo "Cromwell was up but failed to return its version, maybe something went wrong? Continuing waiting..."
      RESULT=1
    fi
  done

  if [ "${RESULT}" -gt "0" ]
  then
    echo "Cromwell never came up after ${ATTEMPTS} attempts"
    exit 1
  else
    echo "Cromwell came up after ${ATTEMPTS} attempts"
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
        clean_up_instance
    fi
}

clean_up_instance() {
    gcloud compute instances delete $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=us-central1-c -q
}

run_test() {
    cd ${CROMWELL_ROOT}

    if [ -z ${GCS_REPORT_PATH} ]
    then
        export GCS_REPORT_PATH="${TEST_CASE_DIRECTORY}/${CROMWELL_VERSION}/${BUILD_ID}"
    fi

    sbt \
      -Dconfig.file=${CENTAUR_CONF_DIR}/centaur.conf \
      -Dcentaur.cromwell.mode="url" \
      -Dcentaur.cromwell.url="http://${CROMWELL_UNDER_TEST}:8000" \
      "centaur/it:testOnly centaur.ExternalTestCaseSpec" | tee centaur.log
    export TEST_RC=$?
}

shutdown() {
    export_logs
    docker-compose -f ${PERF_ROOT}/vm_scripts/docker-compose.yml down
    clean_up
}

read_path_from_vault_json() {
  local path=$1
  local field=$2

  docker run --rm -e VAULT_TOKEN=$VAULT_TOKEN broadinstitute/dsde-toolbox vault read -format=json "${path}" | jq --exit-status -r "${field}"
}


read_service_account_from_vault() {
  mkdir -p mnt
  read_path_from_vault_json "secret/dsp/cromwell/perf/service-account-deployer" '.data.service_account' > mnt/sa.json
}

gcloud_run_as_service_account() {
  local name="$1"
  local command="$2"
  docker run --name $name -v "$(pwd)"/mnt:${DOCKER_ETC_PATH} --rm google/cloud-sdk:slim /bin/bash -c "\
    gcloud auth activate-service-account --key-file ${DOCKER_ETC_PATH}/sa.json 2> /dev/null &&\
    ${command}"
}
