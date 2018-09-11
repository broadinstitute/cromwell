#! bin/bash

echo "Running stress_call_caching_script"

echo "Branch" ${CROMWELL_BRANCH}

# Download workflow files
gsutil cp gs://cromwell-perf-test/workflow_related_files/call_caching_sa2_workflow_options.json workflow_options.json

curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/ss_perf_wf_deploy/scripts/perf/perf_workflow_test/call_caching_stress_test/callCachingStress.wdl -o workflow.wdl
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/ss_perf_wf_deploy/scripts/perf/perf_workflow_test/call_caching_stress_test/callCachingStress_inputs.json -o workflow_inputs.json

curl -X POST "http://localhost:8000/api/workflows/v1" -H "accept: application/json" -H "Content-Type: multipart/form-data" -F "workflowSource=@workflow.wdl" -F "workflowInputs=@workflow_inputs.json;type=application/json" -F "workflowOptions=@workflow_options.json;type=application/json"
