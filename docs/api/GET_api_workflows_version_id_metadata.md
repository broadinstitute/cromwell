This endpoint returns a superset of the data from #get-workflowsversionidlogs in essentially the same format
(i.e. shards are accounted for by an array of maps, in the same order as the shards). 
In addition to shards, every attempt that was made for this call will have its own object as well, in the same order as the attempts.
Workflow metadata includes submission, start, and end datetimes, as well as status, inputs and outputs.
Call-level metadata includes inputs, outputs, start and end datetime, backend-specific job id,
return code, stdout and stderr.  Date formats are ISO with milliseconds.

Accepted parameters are:

* `includeKey` Optional repeated string value, specifies what metadata
  keys to include in the output, matched as a prefix string. Keys that
  are not specified are filtered out. The call keys `attempt` and
  `shardIndex` will always be included. May not be used with
  `excludeKey`.

* `excludeKey` Optional repeated string value, specifies what metadata
  keys to exclude from the output, matched as a prefix string. Keys that
  are specified are filtered out. The call keys `attempt` and
  `shardIndex` will always be included. May not be used with
  `includeKey`.

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata
```

Response:
```
HTTP/1.1 200 OK
Server spray-can/1.3.3 is not blacklisted
Server: spray-can/1.3.3
Date: Thu, 01 Oct 2015 22:18:07 GMT
Content-Type: application/json; charset=UTF-8
Content-Length: 7286
{
  "workflowName": "sc_test",
  "submittedFiles": {
      "inputs": "{}",
      "workflow": "task do_prepare {\n    File input_file\n    command {\n        split -l 1 ${input_file} temp_ && ls -1 temp_?? > files.list\n    }\n    output {\n        Array[File] split_files = read_lines(\"files.list\")\n    }\n}\n# count the number of words in the input file, writing the count to an output file overkill in this case, but simulates a real scatter-gather that would just return an Int (map)\ntask do_scatter {\n    File input_file\n    command {\n        wc -w ${input_file} > output.txt\n    }\n    output {\n        File count_file = \"output.txt\"\n    }\n}\n# aggregate the results back together (reduce)\ntask do_gather {\n    Array[File] input_files\n    command <<<\n        cat ${sep = ' ' input_files} | awk '{s+=$$1} END {print s}'\n    >>>\n    output {\n        Int sum = read_int(stdout())\n    }\n}\nworkflow sc_test {\n    call do_prepare\n    scatter(f in do_prepare.split_files) {\n        call do_scatter {\n            input: input_file = f\n        }\n    }\n    call do_gather {\n        input: input_files = do_scatter.count_file\n    }\n}",
      "options": "{\n\n}",
      "workflowType": "WDL"
  },
  "calls": {
    "sc_test.do_prepare": [
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/stdout",
        "backendStatus": "Done",
        "shardIndex": -1,
        "outputs": {
          "split_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
          ]
        },
        "inputs": {
          "input_file": "/home/jdoe/cromwell/11.txt"
        },
        "runtimeAttributes": {
            "failOnStderr": "true",
            "continueOnReturnCode": "0"
        },
        "callCaching": {
            "allowResultReuse": true,
            "hit": false,
            "result": "Cache Miss",
            "hashes": {
              "output count": "C4CA4238A0B923820DCC509A6F75849B",
              "runtime attribute": {
                "docker": "N/A",
                "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
                "failOnStderr": "68934A3E9455FA72420237EB05902327"
              },
              "output expression": {
                "Array": "D856082E6599CF6EC9F7F42013A2EC4C"
              },
              "input count": "C4CA4238A0B923820DCC509A6F75849B",
              "backend name": "509820290D57F333403F490DDE7316F4",
              "command template": "9F5F1F24810FACDF917906BA4EBA807D",
              "input": {
                "File input_file": "11fa6d7ed15b42f2f73a455bf5864b49"
              }
            },
            "effectiveCallCachingMode": "ReadAndWriteCache"
        },
        "jobId": "34479",
        "returnCode": 0,
        "backend": "Local",
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:55.000-05:00"
      }
    ],
    "sc_test.do_scatter": [
      {
        "executionStatus": "Preempted",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/stdout",
        "backendStatus": "Preempted",
        "shardIndex": 0,
        "outputs": {},
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "callCaching": {
          "allowResultReuse": true,
          "hit": false,
          "result": "Cache Miss",
          "hashes": {
            "output count": "C4CA4238A0B923820DCC509A6F75849B",
            "runtime attribute": {
              "docker": "N/A",
              "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
              "failOnStderr": "68934A3E9455FA72420237EB05902327"
            },
            "output expression": {
              "File count_file": "EF1B47FFA9990E8D058D177073939DF7"
            },
            "input count": "C4CA4238A0B923820DCC509A6F75849B",
            "backend name": "509820290D57F333403F490DDE7316F4",
            "command template": "FD00A1B0AB6A0C97B0737C83F179DDE7",
            "input": {
              "File input_file": "a53794d214dc5dedbcecdf827bf683a2"
            }
          },
         "effectiveCallCachingMode": "ReadAndWriteCache"
        },
        "inputs": {
          "input_file": "f"
        },
        "jobId": "34496",
        "backend": "Local",
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      },
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/stdout",
        "backendStatus": "Done",
        "shardIndex": 0,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt"
        },
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "callCaching": {
          "allowResultReuse": true,
          "hit": false,
          "result": "Cache Miss",
          "hashes": {
            "output count": "C4CA4238A0B923820DCC509A6F75849B",
            "runtime attribute": {
              "docker": "N/A",
              "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
              "failOnStderr": "68934A3E9455FA72420237EB05902327"
            },
            "output expression": {
              "File count_file": "EF1B47FFA9990E8D058D177073939DF7"
            },
            "input count": "C4CA4238A0B923820DCC509A6F75849B",
            "backend name": "509820290D57F333403F490DDE7316F4",
            "command template": "FD00A1B0AB6A0C97B0737C83F179DDE7",
            "input": {
              "File input_file": "a53794d214dc5dedbcecdf827bf683a2"
            }
          },
         "effectiveCallCachingMode": "ReadAndWriteCache"
        },
        "inputs": {
          "input_file": "f"
        },
        "returnCode": 0,
        "jobId": "34965",
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/stderr",
        "attempt": 2,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      },
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/stdout",
        "backendStatus": "Done",
        "shardIndex": 1,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
        },
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "callCaching": {
          "allowResultReuse": true,
          "hit": false,
          "result": "Cache Miss",
          "hashes": {
            "output count": "C4CA4238A0B923820DCC509A6F75849B",
            "runtime attribute": {
              "docker": "N/A",
              "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
              "failOnStderr": "68934A3E9455FA72420237EB05902327"
            },
            "output expression": {
              "File count_file": "EF1B47FFA9990E8D058D177073939DF7"
            },
            "input count": "C4CA4238A0B923820DCC509A6F75849B",
            "backend name": "509820290D57F333403F490DDE7316F4",
            "command template": "FD00A1B0AB6A0C97B0737C83F179DDE7",
            "input": {
              "File input_file": "d3410ade53df34c78488544285cf743c"
            }
          },
         "effectiveCallCachingMode": "ReadAndWriteCache"
        },
        "inputs": {
          "input_file": "f"
        },
        "returnCode": 0,
        "jobId": "34495",
        "backend": "Local",
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      }
    ],
    "sc_test.do_gather": [
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_gather/stdout",
        "backendStatus": "Done",
        "shardIndex": -1,
        "outputs": {
          "sum": 12
        },
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "callCaching": {
          "allowResultReuse": true,
          "hit": false,
          "result": "Cache Miss",
          "hashes": {
            "output count": "C4CA4238A0B923820DCC509A6F75849B",
            "runtime attribute": {
              "docker": "N/A",
              "continueOnReturnCode": "CFCD208495D565EF66E7DFF9F98764DA",
              "failOnStderr": "68934A3E9455FA72420237EB05902327"
            },
            "output expression": {
              "File count_file": "EF1B47FFA9990E8D058D177073939DF7"
            },
            "input count": "C4CA4238A0B923820DCC509A6F75849B",
            "backend name": "509820290D57F333403F490DDE7316F4",
            "command template": "FD00A1B0AB6A0C97B0737C83F179DDE7",
            "input": {
              "File input_file": "e0ef752ab4824939d7947f6012b7c141"
            }
          },
          "effectiveCallCachingMode": "ReadAndWriteCache"
        },
        "inputs": {
          "input_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
          ]
        },
        "returnCode": 0,
        "jobId": "34494",
        "backend": "Local",
        "end": "2016-02-04T13:47:57.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_gather/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      }
    ]
  },
  "outputs": {
    "sc_test.do_gather.sum": 12,
    "sc_test.do_prepare.split_files": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
    ],
    "sc_test.do_scatter.count_file": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
    ]
  },
  "id": "8e592ed8-ebe5-4be0-8dcb-4073a41fe180",
  "inputs": {
    "sc_test.do_prepare.input_file": "/home/jdoe/cromwell/11.txt"
  },
  "labels": {
    "cromwell-workflow-name": "sc_test",
    "cromwell-workflow-id": "cromwell-17633f21-11a9-414f-a95b-2e21431bd67d"
  },
  "submission": "2016-02-04T13:47:55.000-05:00",
  "status": "Succeeded",
  "end": "2016-02-04T13:47:57.000-05:00",
  "start": "2016-02-04T13:47:55.000-05:00"
}
```

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata?includeKey=inputs&includeKey=outputs"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata?includeKey=inputs&includeKey=outputs"
```

Response:
```
HTTP/1.1 200 OK
Server spray-can/1.3.3 is not blacklisted
Server: spray-can/1.3.3
Date: Thu, 01 Oct 2015 22:19:07 GMT
Content-Type: application/json; charset=UTF-8
Content-Length: 4286
{
  "calls": {
    "sc_test.do_prepare": [
      {
        "shardIndex": -1,
        "outputs": {
          "split_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
          ]
        },
        "inputs": {
          "input_file": "/home/jdoe/cromwell/11.txt"
        },
        "attempt": 1
      }
    ],
    "sc_test.do_scatter": [
      {
        "shardIndex": 0,
        "outputs": {},
        "inputs": {
          "input_file": "f"
        },
        "attempt": 1
      },
      {
        "shardIndex": 0,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt"
        },
        "inputs": {
          "input_file": "f"
        },
        "attempt": 2
      },
      {
        "shardIndex": 1,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
        },
        "inputs": {
          "input_file": "f"
        },
        "attempt": 1
      }
    ],
    "sc_test.do_gather": [
      {
        "shardIndex": -1,
        "outputs": {
          "sum": 12
        }
        "inputs": {
          "input_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
          ]
        },
        "attempt": 1
      }
    ]
  },
  "outputs": {
    "sc_test.do_gather.sum": 12,
    "sc_test.do_prepare.split_files": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
    ],
    "sc_test.do_scatter.count_file": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
    ]
  },
  "id": "8e592ed8-ebe5-4be0-8dcb-4073a41fe180",
  "inputs": {
    "sc_test.do_prepare.input_file": "/home/jdoe/cromwell/11.txt"
  }
}
```

The `call` and `workflow` may optionally contain failures shaped like this:
```
"failures": [
  {
    "failure": "The failure message",
    "timestamp": "2016-02-25T10:49:02.066-05:00"
  }
]
```

### Compressing the metadata response

The response from the metadata endpoint can be quite large depending on the workflow. To help with this Cromwell supports gzip encoding the metadata prior to sending it back to the client. In order to enable this, make sure your client is sending the `Accept-Encoding: gzip` header.

For instance, with cURL:
                   
```
$ curl -H "Accept-Encoding: gzip" http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata
```