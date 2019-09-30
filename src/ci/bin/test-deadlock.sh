#!/usr/bin/env bash

set -euo pipefail

pip install docker-py

cromwell_hosts_file="${CROMWELL_BUILD_RESOURCES_DIRECTORY}/cromwell_hosts.txt"

# Get a list of cromwell hosts.
python "${CROMWELL_BUILD_SCRIPTS_DIRECTORY}"/get_cromwell_hosts.py > "${cromwell_hosts_file}"
cromwell_hosts=()
while IFS=$'\n' read -r line; do
    cromwell_hosts+=("$line")
done <"${cromwell_hosts_file}"

# Wait for all of the cromwell hosts to be up.
for host in "${cromwell_hosts[@]}"; do
    "${CROMWELL_BUILD_WAIT_FOR_IT_SCRIPT}" "${host}:8000" -t 180
done

# Submit batches of workflows to the load balancer.
BATCHES_TO_RUN=100
WF_PER_BATCH=10
TOTAL_WFS=$(( BATCHES_TO_RUN * WF_PER_BATCH ))
inputs=$(for ((i=1; i <= $WF_PER_BATCH; i++)); do echo -n ',{}'; done | cut -c2-)
for _ in $(seq $BATCHES_TO_RUN); do
    date
        curl \
        --silent \
        --show-error \
        --request POST \
        --header 'Content-Type: multipart/form-data' \
        --form 'workflowSource=workflow w {call t} task t { command{echo "hello"} }' \
        --form "workflowInputs=[$inputs]" \
        "http://load-balancer:80/api/workflows/v1/batch"
    echo
done

echo "Sleeping up to ten minutes to let workflows run..."
n=0
until [ $n -ge 20 ]
do
    sleep 30
    success=$( \
        curl \
            --silent \
            --show-error \
            --request GET \
            --header "accept: application/json" \
            "http://load-balancer:80/api/workflows/v1/query?status=Succeeded" \
        | grep -Po "totalResultsCount\":([0-9]*)" \
        | cut -d":" -f2 \
        )
    echo "Found $success successful workflows so far"
    [ $success -ge $TOTAL_WFS ] && break
    n=$[$n+1]
done

# Check if all of the hosts are still up.
for host in "${cromwell_hosts[@]}"; do
    echo "Confirming ${host}:8000 is still up"
    curl --fail --silent --show-error "${host}:8000" >/dev/null
done
