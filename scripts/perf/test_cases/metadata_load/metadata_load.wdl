version 1.0

task curl_metadata {
    input {
        String workflowId
        String endpointVersion
    }

    command {
        curl -s 600 http://localhost:8000/api/workflows/~{endpointVersion}/~{workflowId}/metadata > /dev/null
        curl -s 600 http://localhost:8000/api/workflows/~{endpointVersion}/~{workflowId}/metadata > /dev/null
        curl -s 600 http://localhost:8000/api/workflows/~{endpointVersion}/~{workflowId}/metadata > /dev/null
        curl -s 600 http://localhost:8000/api/workflows/~{endpointVersion}/~{workflowId}/metadata > /dev/null
        curl -s 600 http://localhost:8000/api/workflows/~{endpointVersion}/~{workflowId}/metadata > /dev/null
        # Possibly the curls above will timeout, but that's not what we're testing for. We want to check that Cromwell
        # doesn't fall over even if it can't build the metadata under the timeout limit. So exit 0 to make the task succeed 
        # exit 0
    }
}

workflow metadata_load {
    input {
        String endpointVersion
    }
    Array[String] workflowIds = [
        "00b3e3ea-8a9f-4b88-a530-16a155a26db2",
        "02a26bbe-f898-432a-a8c1-4525e54c5389",
        "0315d2d2-3d23-4c7f-a12b-621a1bf850bd",
        "0575401b-7312-4633-981c-453283517050",
        "064169d8-3376-4a5c-9de8-9c8b9ed2e97d",
        "088233c1-262d-470c-bfed-a01d46d56cd5",
        "088f480d-1a8e-4236-997b-e859c6f05775",
        "0a2e7cc0-f871-4dfc-b36c-f888ccabb3dc",
        "0a6592a3-ec33-4c70-8dd8-cf51c2de1379",
        "0c3e5469-669d-4121-9b32-83a6d567416b",
        "0c7a2e57-163a-406d-a9d3-2d6344252c6d",
        "0eec8fa6-840c-47d9-bcfd-bc36de59f90b",
        "0ef72cbe-8154-462d-9de4-7fb9e5c6b5a1",
        "0fcff475-f062-471a-9e09-fe02332c068c",
        "118a62c9-19c0-4dab-a4e6-328dd26447cf"
#        "120c4859-0bf7-47de-b426-01c502ebe6ad",
#        "128ad78d-27b9-4044-a5c4-a8310a628f11",
#        "149edef4-3af1-4773-9373-a7fa4f8b5930",
#        "14f97a3b-d3ac-4ec5-9d85-5ece7d9f7fa7",
#        "15602fa5-d337-4a72-a7e0-ec3766fc0450",
#        "17262d7a-addc-45c2-bd02-ad0d612a802d",
#        "17a9f931-5194-4948-909b-ee1c5b8449b3",
#        "17b76b73-682e-46c6-871f-dacf0fb6a99b",
#        "1b0573e5-719f-4566-969b-c271528eecdc",
#        "1be5d749-765f-4bd6-8605-c0bfce2ae4b4",
#        "21b297f6-643a-4c75-8000-38ae57ef642f",
#        "23766f7c-1459-45d6-808b-c6c7343b680c",
#        "2551d04c-38ed-424b-a614-fa9900aaf808",
#        "29da3bf9-959b-4308-b504-63c2d3794a98",
#        "2b75c685-76ff-49a0-802a-ef0eddc4b793"
    ]
    
    scatter(id in workflowIds) {
        call curl_metadata { input: workflowId = id, endpointVersion = endpointVersion }
    }
}
