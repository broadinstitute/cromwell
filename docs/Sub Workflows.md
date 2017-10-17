WDL allows the execution of an entire workflow as a step in a larger workflow (see WDL SPEC for more details), which is what will be referred to as a sub workflow going forward.
Cromwell supports execution of such workflows. Note that sub workflows can themselves contain sub workflows, etc... There is no limitation as to how deeply workflows can be nested.

## Execution

Sub workflows are executed exactly as a task would be.
*This means that if another call depends on an output of a sub workflow, this call will run when the whole sub workflow completes (successfully).*
For example, in the following case :

`main.wdl`
```
import "sub_wdl.wdl" as sub

workflow main_workflow {

    call sub.hello_and_goodbye { input: hello_and_goodbye_input = "sub world" }
    
    # call myTask { input: hello_and_goodbye.hello_output }
    
    output {
        String main_output = hello_and_goodbye.hello_output
    }
}
```

`sub_wdl.wdl`
```
task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  runtime {
      docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string(stdout())
  }
}

task goodbye {
  String addressee
  command {
    echo "Goodbye ${addressee}!"
  }
  runtime {
      docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string(stdout())
  }
}

workflow hello_and_goodbye {
  String hello_and_goodbye_input
  
  call hello {input: addressee = hello_and_goodbye_input }
  call goodbye {input: addressee = hello_and_goodbye_input }
  
  output {
    String hello_output = hello.salutation
    String goodbye_output = goodbye.salutation
  }
}
```

`myTask` will start only when hello_and_goodbye completes (which means all of its calls are done), even though `myTask` only needs the output of hello in the hello_and_goodbye sub workflow. 
If hello_and_goodbye fails, then `myTask` won't be executed.
Only workflow outputs are visible outside a workflow, which means that references to outputs produced by a sub workflow will only be valid if those outputs are exposed in the workflow output section.

Sub workflows are executed in the context of a main workflow, which means that operations that are normally executed once per workflow (set up, clean up, outputs copying, log copying, etc...)
will NOT be re-executed for each sub workflow. For instance if a resource is created during workflow initialization, sub workflows will need to share this same resource.
Workflow outputs will be copied for the main root workflow but not for intermediate sub workflows.

Restarts, aborts, and call-caching work exactly as they would with tasks. 
All tasks run by a sub workflow are eligible for call caching under the same rules as any other task.
However, workflows themselves are not cached as such. Which means that running the exact same workflow twice with call caching on will trigger each task to cache individually,
but not the workflow itself.

The root path for sub workflow execution files (scripts, output files, logs) will be under the parent workflow call directory.
For example, the execution directory for the above main workflow would look like the following:

```
cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/ <- main workflow id
└── call-hello_and_goodbye <- call directory for call hello_and_goodbye in the main workflow
    └── hello_and_goodbye <- name of the sub workflow 
        └── a6365f91-c807-465a-9186-a5d3da98fe11 <- sub workflow id
            ├── call-goodbye
            │   └── execution
            │       ├── rc
            │       ├── script
            │       ├── script.background
            │       ├── script.submit
            │       ├── stderr
            │       ├── stderr.background
            │       ├── stdout
            │       └── stdout.background
            └── call-hello
                └── execution
                    ├── rc
                    ├── script
                    ├── script.background
                    ├── script.submit
                    ├── stderr
                    ├── stderr.background
                    ├── stdout
                    └── stdout.background

```

## Metadata
Each sub workflow will have its own workflow ID. This ID will appear in the metadata of the parent workflow, in the call section corresponding to the sub workflow, under the "subWorkflowId" attribute.
For example, querying the `main_workflow` metadata above (minus the `myTask` call) , could result in something like this:

`GET /api/workflows/v2/1d919bd4-d046-43b0-9918-9964509689dd/metadata`

```
{
  "workflowName": "main_workflow",
  "submittedFiles": {
    "inputs": "{}",
    "workflow": "import \"sub_wdl.wdl\" as sub\n\nworkflow main_workflow {\n\n    call sub.hello_and_goodbye { input: hello_and_goodbye_input = \"sub world\" }\n    \n    # call myTask { input: hello_and_goodbye.hello_output }\n    \n    output {\n        String main_output = hello_and_goodbye.hello_output\n    }\n}",
    "options": "{\n\n}"
  },
  "calls": {
    "main_workflow.hello_and_goodbye": [
      {
        "executionStatus": "Done",
        "shardIndex": -1,
        "outputs": {
          "goodbye_output": "Goodbye sub world!",
          "hello_output": "Hello sub world!"
        },
        "inputs": {
          "hello_and_goodbye_input": "sub world"
        },
        "end": "2016-11-17T14:13:41.117-05:00",
        "attempt": 1,
        "start": "2016-11-17T14:13:39.236-05:00",
        "subWorkflowId": "a6365f91-c807-465a-9186-a5d3da98fe11"
      }
    ]
  },
  "outputs": {
    "main_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd",
  "id": "1d919bd4-d046-43b0-9918-9964509689dd",
  "inputs": {},
  "submission": "2016-11-17T14:13:39.104-05:00",
  "status": "Succeeded",
  "end": "2016-11-17T14:13:41.120-05:00",
  "start": "2016-11-17T14:13:39.204-05:00"
}
```

The sub workflow ID can be queried separately:

`GET /api/workflows/v2/a6365f91-c807-465a-9186-a5d3da98fe11/metadata`

```
{
  "workflowName": "hello_and_goodbye",
  "calls": {
    "sub.hello_and_goodbye.hello": [
      {
        "executionStatus": "Done",
        "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stdout",
        "shardIndex": -1,
        "outputs": {
          "salutation": "Hello sub world!"
        },
        "runtimeAttributes": {
          "docker": "ubuntu:latest",
          "failOnStderr": false,
          "continueOnReturnCode": "0"
        },
        "cache": {
          "allowResultReuse": true
        },
        "Effective call caching mode": "CallCachingOff",
        "inputs": {
          "addressee": "sub world"
        },
        "returnCode": 0,
        "jobId": "49830",
        "backend": "Local",
        "end": "2016-11-17T14:13:40.712-05:00",
        "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stderr",
        "callRoot": "/cromwell/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello",
        "attempt": 1,
        "executionEvents": [
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "Pending",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "RequestingExecutionToken",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "PreparingJob",
            "endTime": "2016-11-17T14:13:39.243-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.243-05:00",
            "description": "RunningJob",
            "endTime": "2016-11-17T14:13:40.704-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:40.704-05:00",
            "description": "UpdatingJobStore",
            "endTime": "2016-11-17T14:13:40.712-05:00"
          }
        ],
        "start": "2016-11-17T14:13:39.239-05:00"
      }
    ],
    "sub.hello_and_goodbye.goodbye": [
      {
        "executionStatus": "Done",
        "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stdout",
        "shardIndex": -1,
        "outputs": {
          "salutation": "Goodbye sub world!"
        },
        "runtimeAttributes": {
          "docker": "ubuntu:latest",
          "failOnStderr": false,
          "continueOnReturnCode": "0"
        },
        "cache": {
          "allowResultReuse": true
        },
        "Effective call caching mode": "CallCachingOff",
        "inputs": {
          "addressee": "sub world"
        },
        "returnCode": 0,
        "jobId": "49831",
        "backend": "Local",
        "end": "2016-11-17T14:13:41.115-05:00",
        "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stderr",
        "callRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye",
        "attempt": 1,
        "executionEvents": [
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "Pending",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "RequestingExecutionToken",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "PreparingJob",
            "endTime": "2016-11-17T14:13:39.243-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.243-05:00",
            "description": "RunningJob",
            "endTime": "2016-11-17T14:13:41.112-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:41.112-05:00",
            "description": "UpdatingJobStore",
            "endTime": "2016-11-17T14:13:41.115-05:00"
          }
        ],
        "start": "2016-11-17T14:13:39.239-05:00"
      }
    ]
  },
  "outputs": {
    "goodbye_output": "Goodbye sub world!",
    "hello_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11",
  "id": "a6365f91-c807-465a-9186-a5d3da98fe11",
  "inputs": {
    "hello_and_goodbye_input": "sub world"
  },
  "status": "Succeeded",
  "parentWorkflowId": "1d919bd4-d046-43b0-9918-9964509689dd",
  "end": "2016-11-17T14:13:41.116-05:00",
  "start": "2016-11-17T14:13:39.236-05:00"
}
```

It's also possible to set the URL query parameter `expandSubWorkflows` to `true` to automatically include sub workflows metadata (`false` by default).

`GET api/workflows/v2/1d919bd4-d046-43b0-9918-9964509689dd/metadata?expandSubWorkflows=true`

```
{
  "workflowName": "main_workflow",
  "submittedFiles": {
    "inputs": "{}",
    "workflow": "import \"sub_wdl.wdl\" as sub\n\nworkflow main_workflow {\n\n    call sub.hello_and_goodbye { input: hello_and_goodbye_input = \"sub world\" }\n    \n    # call myTask { input: hello_and_goodbye.hello_output }\n    \n    output {\n        String main_output = hello_and_goodbye.hello_output\n    }\n}",
    "options": "{\n\n}"
  },
  "calls": {
    "main_workflow.hello_and_goodbye": [{
      "executionStatus": "Done",
      "subWorkflowMetadata": {
        "workflowName": "hello_and_goodbye",
        "calls": {
          "sub.hello_and_goodbye.hello": [{
            "executionStatus": "Done",
            "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stdout",
            "shardIndex": -1,
            "outputs": {
              "salutation": "Hello sub world!"
            },
            "runtimeAttributes": {
              "docker": "ubuntu:latest",
              "failOnStderr": false,
              "continueOnReturnCode": "0"
            },
            "cache": {
              "allowResultReuse": true
            },
            "Effective call caching mode": "CallCachingOff",
            "inputs": {
              "addressee": "sub world"
            },
            "returnCode": 0,
            "jobId": "49830",
            "backend": "Local",
            "end": "2016-11-17T14:13:40.712-05:00",
            "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stderr",
            "callRoot": "cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello",
            "attempt": 1,
            "executionEvents": [{
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "Pending",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "RequestingExecutionToken",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "PreparingJob",
              "endTime": "2016-11-17T14:13:39.243-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.243-05:00",
              "description": "RunningJob",
              "endTime": "2016-11-17T14:13:40.704-05:00"
            }, {
              "startTime": "2016-11-17T14:13:40.704-05:00",
              "description": "UpdatingJobStore",
              "endTime": "2016-11-17T14:13:40.712-05:00"
            }],
            "start": "2016-11-17T14:13:39.239-05:00"
          }],
          "sub.hello_and_goodbye.goodbye": [{
            "executionStatus": "Done",
            "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stdout",
            "shardIndex": -1,
            "outputs": {
              "salutation": "Goodbye sub world!"
            },
            "runtimeAttributes": {
              "docker": "ubuntu:latest",
              "failOnStderr": false,
              "continueOnReturnCode": "0"
            },
            "cache": {
              "allowResultReuse": true
            },
            "Effective call caching mode": "CallCachingOff",
            "inputs": {
              "addressee": "sub world"
            },
            "returnCode": 0,
            "jobId": "49831",
            "backend": "Local",
            "end": "2016-11-17T14:13:41.115-05:00",
            "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stderr",
            "callRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye",
            "attempt": 1,
            "executionEvents": [{
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "Pending",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "RequestingExecutionToken",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "PreparingJob",
              "endTime": "2016-11-17T14:13:39.243-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.243-05:00",
              "description": "RunningJob",
              "endTime": "2016-11-17T14:13:41.112-05:00"
            }, {
              "startTime": "2016-11-17T14:13:41.112-05:00",
              "description": "UpdatingJobStore",
              "endTime": "2016-11-17T14:13:41.115-05:00"
            }],
            "start": "2016-11-17T14:13:39.239-05:00"
          }]
        },
        "outputs": {
          "goodbye_output": "Goodbye sub world!",
          "hello_output": "Hello sub world!"
        },
        "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11",
        "id": "a6365f91-c807-465a-9186-a5d3da98fe11",
        "inputs": {
          "hello_and_goodbye_input": "sub world"
        },
        "status": "Succeeded",
        "parentWorkflowId": "1d919bd4-d046-43b0-9918-9964509689dd",
        "end": "2016-11-17T14:13:41.116-05:00",
        "start": "2016-11-17T14:13:39.236-05:00"
      },
      "shardIndex": -1,
      "outputs": {
        "goodbye_output": "Goodbye sub world!",
        "hello_output": "Hello sub world!"
      },
      "inputs": {
        "hello_and_goodbye_input": "sub world"
      },
      "end": "2016-11-17T14:13:41.117-05:00",
      "attempt": 1,
      "start": "2016-11-17T14:13:39.236-05:00"
    }]
  },
  "outputs": {
    "main_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd",
  "id": "1d919bd4-d046-43b0-9918-9964509689dd",
  "inputs": {

  },
  "submission": "2016-11-17T14:13:39.104-05:00",
  "status": "Succeeded",
  "end": "2016-11-17T14:13:41.120-05:00",
  "start": "2016-11-17T14:13:39.204-05:00"
}
```