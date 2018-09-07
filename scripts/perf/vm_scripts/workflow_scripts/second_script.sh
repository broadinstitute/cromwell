#! bin/bash

echo "Inside second script"

mkdir /second_script_folder
cd /second_script_folder

curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/develop/centaur/src/main/resources/standardTestCases/hello/hello.wdl -o workflow.wdl
curl -L https://raw.githubusercontent.com/broadinstitute/cromwell/develop/centaur/src/main/resources/standardTestCases/hello/hello.inputs -o workflow_inputs.json

curl -X POST "http://localhost:8000/api/workflows/v1" -H "accept: application/json" -H "Content-Type: multipart/form-data" -F "workflowSource=@workflow.wdl" -F "workflowInputs=@workflow_inputs.json;type=application/json"
