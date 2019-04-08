#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

source "${BASH_SOURCE%/*}/test.inc.sh" || source test.inc.sh

# A set of functions for use in the Kubernetes test script.

GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON="cromwell-centaur-service-account.json"
GOOGLE_ZONE=us-central1-c
DOCKER_ETC_PATH=/usr/share/etc

cromwell::build::setup_common_environment

cromwell::build::setup_centaur_environment

GOOGLE_PROJECT=$(cat "$CROMWELL_BUILD_RESOURCES_DIRECTORY/$GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON" | jq -r .project_id)

# Takes a single string argument and `echo`s a possibly modified version of that argument with non-alphanumeric
# characters converted to dashes. TODO: restrict the initial character as necessary
cromwell::kube::google_safe_name() {
  echo -n "$1" | tr -c '[[:digit:][:alpha:]]' '-'
}

cromwell::kube::google_safe_build_name() {
  echo -n "$(cromwell::kube::google_safe_name ${CROMWELL_BUILD_PROVIDER}-${CROMWELL_BUILD_NUMBER:-$RANDOM})"
}

# Creates a Google friendly identifier name specific to this build based on a single argument.
cromwell::kube::centaur_gke_name() {
  local prefix="centaur-gke"
  local build_name="$(cromwell::kube::google_safe_build_name)"
  local arg=$1
  echo -n "${prefix}-${arg}-${build_name}"
}

# Run a specified command after activating the specified service account.
#
# Usage: cromwell::kube::gcloud_run_as_service_account command
cromwell::kube::gcloud_run_as_service_account() {
  local command="$1"
  docker run -v "$CROMWELL_BUILD_RESOURCES_DIRECTORY:$DOCKER_ETC_PATH" -e DOCKER_ETC_PATH --rm google/cloud-sdk:latest /bin/bash -c "\
    gcloud auth activate-service-account --key-file $DOCKER_ETC_PATH/${GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON} 2> /dev/null && $command "
}

# Configure kubectl for gke then run the specified command
cromwell::kube::gcloud_run_kubectl_command_as_service_account() {
  local gkeClusterName="$1"
  local command="$2"
  cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT container clusters get-credentials --zone $GOOGLE_ZONE $gkeClusterName && $command"
}

cromwell::kube::generate_cloud_sql_instance_name() {
  echo -n $(cromwell::kube::centaur_gke_name "cloudsql")
}

# Starts a Cloud SQL instance with the specified name.
cromwell::kube::create_cloud_sql_instance() {
  local cloudSqlInstanceName="$1"
  local cloudSqlPassword="$(cat ${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell-centaur-gke-cloudsql.json | jq -r '.db_pass')"

  # Kick off async creation of the Cloud SQL instance (synchronous creation fails sometimes when the cloud is angry).
  local operationName=$(cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT sql instances create --async --zone $GOOGLE_ZONE --storage-size=10GB --database-version=MYSQL_5_7 $cloudSqlInstanceName --format='value(name)'")

  # wait for it
  cromwell::kube::gcloud_run_as_service_account \
    "gcloud beta sql operations wait --timeout=900 --project $GOOGLE_PROJECT ${operationName}"

  # Create a user.
  cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT sql users create cromwell --instance $cloudSqlInstanceName --password='${cloudSqlPassword}'"

  # Create the cromwell test database.
  cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT sql databases create cromwell_test --instance $cloudSqlInstanceName"
}

cromwell::kube::destroy_cloud_sql_instance() {
  local instanceName="$1"
  cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT --quiet sql instances delete $instanceName"
}

# Returns the connection name for the specific Cloud SQL instance name.
#
# Usage: cromwell::kube::connection_name_for_cloud_sql_instance instance_name
cromwell::kube::connection_name_for_cloud_sql_instance() {
  # TOL It appears the connectionName can be inferred (<project>:<region>:<instance name>), so it may not be necessary to query.
  local instanceName="$1"
  # TODO it would be nice to make this async like the Cloud SQL creation but the `gcloud container clusters create --async` doesn't
  # TODO seem to return an operation ID. It might be possible to list operations and search for the known cluster name but yuck.
  echo -n $(cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT sql instances describe $instanceName --format='value(connectionName)'" | tr -d '\n')
}

cromwell::kube::generate_gke_cluster_name() {
  echo -n $(cromwell::kube::centaur_gke_name "cluster")
}

# Create a GKE cluster with the specified name.
cromwell::kube::create_gke_cluster() {
  local gkeClusterName="$1"
  cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT container clusters create --zone $GOOGLE_ZONE $gkeClusterName --machine-type=n1-standard-2 --num-nodes=3"

  echo -n ${gkeClusterName}
}

cromwell:kube::destroy_gke_cluster() {
  local gkeClusterName="$1"
    cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT --quiet container clusters delete $gkeClusterName --zone $GOOGLE_ZONE"
}

cromwell::kube::generate_gcr_tag() {
  local buildName=$(cromwell::kube::google_safe_build_name)
  echo -n "gcr.io/$GOOGLE_PROJECT/centaur-gke/cromwell:${buildName}"
}

cromwell::kube::tag_for_gcr() {
  local image="$1"
  local tag="$2"
  docker tag ${image} ${tag}
}

cromwell::kube::gcr_login() {
  cat "$CROMWELL_BUILD_RESOURCES_DIRECTORY/$GOOGLE_CENTAUR_SERVICE_ACCOUNT_JSON" | docker login -u _json_key --password-stdin https://gcr.io
}

cromwell::kube::push_to_gcr() {
  local dummy="$1"
  local tag="$2"
  docker tag broadinstitute/cromwell:${dummy} ${tag}
  docker push ${tag}
}

cromwell::kube::delete_from_gcr() {
  local tag="$1"
  cromwell::kube::gcloud_run_as_service_account \
    "gcloud --project $GOOGLE_PROJECT --quiet container images delete $tag"
}

cromwell::kube::create_secrets() {
  local from_files=""
  for file in ${CROMWELL_BUILD_RESOURCES_DIRECTORY}/*.conf ${CROMWELL_BUILD_RESOURCES_DIRECTORY}/*.json
  do
    # This is going to run inside the gcloud Docker container which mounts the resources directory at $DOCKER_ETC_PATH
    from_files+="--from-file=${DOCKER_ETC_PATH}/$(basename ${file}) "
  done

  # In the final version of this Cromwell secrets will not need to be named in a build-specific way because they will be
  # scoped to build-specific clusters. But for testing the cluster is currently fixed, so this better not exist beforehand.
  local command="kubectl create secret generic cromwell-secrets ${from_files}"
  echo "Creating secrets with command: $command"

  cromwell::kube::gcloud_run_kubectl_command_as_service_account \
    "${KUBE_CLUSTER_NAME}" "${command}"
}

cromwell::kube::private::wait_for_summarizer_deployment() {
  local i=0
  local maxRetries=30
  local waitTime=10

  while [[ "${i}" -le "${maxRetries}" ]]
  do
    replicaCount=$(cromwell::kube::gcloud_run_kubectl_command_as_service_account  "${KUBE_CLUSTER_NAME}" \
      "kubectl get deployment cromwell-summarizer -o jsonpath='{.status.availableReplicas}'")

    if [[ "${replicaCount}" -eq "1" ]]; then
      echo "Summarizer deployment says it's ready, but let's sleep a little more just to be sure"
      sleep 60
      break
    fi

    i=$((i+1))
    echo "Summarizer deployment not ready on iteration ${i}, sleeping ${waitTime} seconds"
    sleep ${waitTime}
  done
}

cromwell::kube::get_load_balancer_ip() {
  local outfile="$1"
  local i=0
  local maxRetries=30
  local waitTime=10
  local ip=""

  while [[ "${i}" -le "${maxRetries}" ]]
  do
    ip=$(cromwell::kube::gcloud_run_kubectl_command_as_service_account  "${KUBE_CLUSTER_NAME}" \
      "kubectl get services cromwell-frontend -o jsonpath='{.status.loadBalancer.ingress[0].ip}' ")

    if [[ -n "${ip}" ]]; then
      echo "Load balancer IP found"
      echo "${ip}" > "${outfile}"
      break
    fi

    i=$((i+1))
    echo "Front end load balancer IP not available on iteration ${i}, sleeping ${waitTime} seconds"
    sleep ${waitTime}
  done
}


cromwell::kube::start_cromwell() {
  for instance_type in summarizer frontend backend
  do
    KUBE_CROMWELL_INSTANCE_TYPE="${instance_type}"
    cromwell::kube::gcloud_run_kubectl_command_as_service_account \
      "${KUBE_CLUSTER_NAME}" "kubectl apply -f ${DOCKER_ETC_PATH}/${instance_type}-cromwell-deployment.yaml"
    if [[ "$instance_type" == "summarizer" ]]; then
      # Have the sole summarizer instance come up first and run Liquibase before applying the other deployments
      cromwell::kube::private::wait_for_summarizer_deployment
    fi
  done
}

cromwell::kube::private::build_environment_variable_arguments() {
  local ret=""
  for e in $(compgen -v | grep '^KUBE_')
  do
    value=$(eval echo "\${${e}}")
    ret+=" --env ${e}=${value} "
  done
  echo ${ret}
}

cromwell::kube::create_deployment_files() {
  local dtmpl="${CROMWELL_BUILD_RESOURCES_SOURCES}/cromwell-deployment.yaml.dtmpl"

  for type in frontend backend summarizer
  do
    KUBE_CROMWELL_INSTANCE_TYPE="${type}"
    if [[ "${type}" == "frontend" ]]; then
      KUBE_DEPLOYMENT_REPLICA_COUNT=2
    elif [[ "${type}" == "backend" ]]; then
      KUBE_DEPLOYMENT_REPLICA_COUNT=3
    elif [[ "${type}" == "summarizer" ]]; then
      KUBE_DEPLOYMENT_REPLICA_COUNT=1
    fi

    local env=$(cromwell::kube::private::build_environment_variable_arguments)

    local outfile="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/${type}-cromwell-deployment.yaml"
    local command="docker run --rm ${env} -v ${CROMWELL_BUILD_RESOURCES_SOURCES}:${CROMWELL_BUILD_RESOURCES_SOURCES} -v ${CROMWELL_BUILD_RESOURCES_DIRECTORY}:${CROMWELL_BUILD_RESOURCES_DIRECTORY} -it broadinstitute/dsde-toolbox consul-template -once -template=${dtmpl}:${outfile}"
    echo "deployment template: ${command}"
    eval ${command}
  done
}
