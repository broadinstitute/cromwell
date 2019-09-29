#!/usr/bin/env bash

set -euo pipefail

pip install docker-py

wget https://raw.githubusercontent.com/vishnubob/wait-for-it/db04971/wait-for-it.sh
chmod +x wait-for-it.sh

# Get a list of cromwell hosts.
python /cromwell/bin/get_cromwell_hosts.py > cromwell_hosts.txt
cromwell_hosts=()
while IFS=$'\n' read -r line; do
  cromwell_hosts+=("$line")
done < cromwell_hosts.txt

# Wait for all of the cromwell hosts to be up.
for host in "${cromwell_hosts[@]}"; do
  ./wait-for-it.sh "${host}:8000" -t 180
done

# Submit batches of workflows to the load balancer.
BATCHES_TO_RUN=100
WF_PER_BATCH=10
TOTAL_WFS=$(( BATCHES_TO_RUN * WF_PER_BATCH ))
inputs=$(for ((i=1; i <= $WF_PER_BATCH; i++)); do echo -n ',{}'; done | cut -c2-)
for _ in $(seq $BATCHES_TO_RUN); do
  echo $(date)
  curl \
    -s -X POST \
    "http://lb:80/api/workflows/v1/batch" \
    -H 'Content-Type: multipart/form-data' \
    -F 'workflowSource=workflow w {call t} task t { command{echo "hello"} }' \
    -F "workflowInputs=[$inputs]"
  echo
done

echo "Sleeping up to ten minutes to let workflows run..."
n=0
until [ $n -ge 20 ]
do
    sleep 30
    success=$(curl -s -X GET "http://lb:80/api/workflows/v1/query?status=Succeeded" -H "accept: application/json" | grep -Po "totalResultsCount\":([0-9]*)" | cut -d":" -f2)
    echo "Found $success successful workflows so far"
    [ $success -ge $TOTAL_WFS ] && break
    n=$[$n+1]
done

# Check if all of the hosts are still up.
for host in "${cromwell_hosts[@]}"; do
  wget "${host}:8000"
done
