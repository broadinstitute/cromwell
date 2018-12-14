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
for _ in $(seq 50); do
  curl \
    -X POST \
    "http://lb:80/api/workflows/v1/batch" \
    -H 'Content-Type: multipart/form-data' \
    -F 'workflowSource=workflow w {call t} task t { command{echo "hello"} }' \
    -F 'workflowInputs=[{},{},{},{},{},{},{},{},{},{}]'
done

echo "Sleeping five minutes to let workflows run..."
sleep 300

# Check if all of the hosts are still up.
for host in "${cromwell_hosts[@]}"; do
  wget "${host}:8000"
done
