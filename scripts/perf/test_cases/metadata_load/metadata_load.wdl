version 1.0

task install_ab {
    command {
        apt-get install --assume-yes apache2-utils
    }
    output {
        Boolean done = true
    }
}

task metadata_ab {
    input {
        Boolean start
        String workflowId
        String endpointVersion
        Int nbRequests
        Int concurrency
        Int timeout = 70
    }

    command {
        ab -n ~{nbRequests} -c ~{concurrency} -s ~{timeout} http://localhost:8000/api/workflows/~{endpointVersion}/~{workflowId}/metadata
    }
}

workflow metadata_load {
    input {
        String endpointVersion
    }

    # 100 000 rows each
    Array[String] largeWorkflowIds = [
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
    ]
    
    # 1500 rows
    String smallWorkflowId = "7a35b251-8083-481a-9089-085e04314258"
    
    call install_ab

    scatter(id in largeWorkflowIds) {
        call metadata_ab { input: workflowId = id, endpointVersion = endpointVersion, start = install_ab.done, nbRequests = 1000, concurrency = 50 }
    }
    
    call metadata_ab { input: workflowId = smallWorkflowId, endpointVersion = endpointVersion, start = install_ab.done, nbRequests = 10000, concurrency = 100 }
}
