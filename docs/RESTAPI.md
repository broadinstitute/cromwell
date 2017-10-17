The `server` subcommand on the executable JAR will start an HTTP server which can accept workflow files to run as well as check status and output of existing workflows.

The following sub-sections define which HTTP Requests the web server can accept and what they will return.  Example HTTP requests are given in [HTTPie](https://github.com/jkbrzt/httpie) and [cURL](https://curl.haxx.se/)

## REST API Versions

All web server requests include an API version in the url. The current version is `v1`.

## POST /api/workflows/:version

This endpoint accepts a POST request with a `multipart/form-data` encoded body.  The form fields that may be included are:

* `workflowSource` - *Required* Contains the workflow source file to submit for execution.
* `workflowType` - *Optional* The type of the `workflowSource`.
   If not specified, returns the `workflow-options.default.workflow-type` configuration value with a default of `WDL`.
* `workflowTypeVersion` - *Optional* The version of the `workflowType`, for example "draft-2".
   If not specified, returns the `workflow-options.default.workflow-type-version` configuration value with no default.
* `workflowInputs` - *Optional* JSON file containing the inputs.  For WDL workflows a skeleton file can be generated from [wdltool](https://github.com/broadinstitute/wdltool) using the "inputs" subcommand.
* `workflowInputs_n` - *Optional* Where `n` is an integer. JSON file containing the 'n'th set of auxiliary inputs.
* `workflowOptions` - *Optional* JSON file containing options for this workflow execution.  See the [run](#run) CLI sub-command for some more information about this.
* `customLabels` - *Optional* JSON file containing a set of custom labels to apply to this workflow. See [Labels](#labels) for the expected format.
* `workflowDependencies` - *Optional* ZIP file containing workflow source files that are used to resolve import statements.

Regarding the workflowInputs parameter, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2.
Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file.

Additionally, although Swagger has a limit of 5 JSON input files, the REST endpoint itself can accept an unlimited number of JSON input files.


cURL:

```
$ curl -v "localhost:8000/api/workflows/v1" -F workflowSource=@src/main/resources/3step.wdl -F workflowInputs=@test.json
```

HTTPie:

```
$ http --print=hbHB --form POST localhost:8000/api/workflows/v1 workflowSource=@src/main/resources/3step.wdl workflowInputs@inputs.json
```

Request:

```
POST /api/workflows/v1 HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 730
Content-Type: multipart/form-data; boundary=64128d499e9e4616adea7d281f695dca
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--64128d499e9e4616adea7d281f695dca
Content-Disposition: form-data; name="workflowSource"

task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}

task cgrep {
  command {
    grep '${pattern}' ${File in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

task wc {
  command {
    cat ${File in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

workflow three_step {
  call ps
  call cgrep {
    input: in_file=ps.procs
  }
  call wc {
    input: in_file=ps.procs
  }
}

--64128d499e9e4616adea7d281f695dca
Content-Disposition: form-data; name="workflowInputs"; filename="inputs.json"

{
    "three_step.cgrep.pattern": "..."
}

--64128d499e9e4616adea7d281f695dca--
```

Response:

```
HTTP/1.1 201 Created
Content-Length: 74
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:28 GMT
Server: spray-can/1.3.3

{
    "id": "69d1d92f-3895-4a7b-880a-82535e9a096e",
    "status": "Submitted"
}
```

To specify workflow options as well:


cURL:

```
$ curl -v "localhost:8000/api/workflows/v1" -F workflowSource=@wdl/jes0.wdl -F workflowInputs=@wdl/jes0.json -F workflowOptions=@options.json
```

HTTPie:

```
http --print=HBhb --form POST http://localhost:8000/api/workflows/v1 workflowSource=@wdl/jes0.wdl workflowInputs@wdl/jes0.json workflowOptions@options.json
```

Request (some parts truncated for brevity):

```
POST /api/workflows/v1 HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 1472
Content-Type: multipart/form-data; boundary=f3fd038395644de596c460257626edd7
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="workflowSource"

task x { ... }
task y { ... }
task z { ... }

workflow myworkflow {
  call x
  call y
  call z {
    input: example="gs://my-bucket/cromwell-executions/myworkflow/example.txt", int=3000
  }
}

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="workflowInputs"; filename="jes0.json"

{
  "myworkflow.x.x": "100"
}

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="workflowOptions"; filename="options.json"

{
  "jes_gcs_root": "gs://myworkflow-dev/workflows"
}

--f3fd038395644de596c460257626edd7--
```

## POST /api/workflows/:version/batch

This endpoint accepts a POST request with a `multipart/form-data`
encoded body.  The form fields that may be included are:

* `workflowSource` - *Required* Contains the workflow source file to submit for
execution.
* `workflowInputs` - *Required* JSON file containing the inputs in a
JSON array. For WDL workflows a skeleton file for a single inputs json element can be
generated from [wdltool](https://github.com/broadinstitute/wdltool)
using the "inputs" subcommand. The orderded endpoint responses will
contain one workflow submission response for each input, respectively.
* `workflowOptions` - *Optional* JSON file containing options for this
workflow execution.  See the [run](#run) CLI sub-command for some more
information about this.
* `workflowDependencies` - *Optional* ZIP file containing workflow source files that are used to resolve import statements. Applied equally to all workflowInput sets.

cURL:

```
$ curl -v "localhost:8000/api/workflows/v1/batch" -F workflowSource=@src/main/resources/3step.wdl -F workflowInputs=@test_array.json
```

HTTPie:

```
$ http --print=hbHB --form POST localhost:8000/api/workflows/v1/batch workflowSource=@src/main/resources/3step.wdl workflowInputs@inputs_array.json
```

Request:

```
POST /api/workflows/v1/batch HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 750
Content-Type: multipart/form-data; boundary=64128d499e9e4616adea7d281f695dcb
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--64128d499e9e4616adea7d281f695dcb
Content-Disposition: form-data; name="workflowSource"

task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}

task cgrep {
  command {
    grep '${pattern}' ${File in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

task wc {
  command {
    cat ${File in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

workflow three_step {
  call ps
  call cgrep {
    input: in_file=ps.procs
  }
  call wc {
    input: in_file=ps.procs
  }
}

--64128d499e9e4616adea7d281f695dcb
Content-Disposition: form-data; name="workflowInputs"; filename="inputs_array.json"

[
    {
        "three_step.cgrep.pattern": "..."
    },
    {
        "three_step.cgrep.pattern": "..."
    }
]

--64128d499e9e4616adea7d281f695dcb--
```

Response:

```
HTTP/1.1 201 Created
Content-Length: 96
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:28 GMT
Server: spray-can/1.3.3

[
    {
        "id": "69d1d92f-3895-4a7b-880a-82535e9a096e",
        "status": "Submitted"
    },
    {
        "id": "69d1d92f-3895-4a7b-880a-82535e9a096f",
        "status": "Submitted"
    }
]
```

To specify workflow options as well:


cURL:

```
$ curl -v "localhost:8000/api/workflows/v1/batch" -F workflowSource=@wdl/jes0.wdl -F workflowInputs=@wdl/jes0_array.json -F workflowOptions=@options.json
```

HTTPie:

```
http --print=HBhb --form POST http://localhost:8000/api/workflows/v1/batch workflowSource=@wdl/jes0.wdl workflowInputs@wdl/jes0_array.json workflowOptions@options.json
```

Request (some parts truncated for brevity):

```
POST /api/workflows/v1/batch HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 1492
Content-Type: multipart/form-data; boundary=f3fd038395644de596c460257626edd8
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--f3fd038395644de596c460257626edd8
Content-Disposition: form-data; name="workflowSource"

task x { ... }
task y { ... }
task z { ... }

workflow myworkflow {
  call x
  call y
  call z {
    input: example="gs://my-bucket/cromwell-executions/myworkflow/example.txt", int=3000
  }
}

--f3fd038395644de596c460257626edd8
Content-Disposition: form-data; name="workflowInputs"; filename="jes0_array.json"

[
  {
    "myworkflow.x.x": "100"
  }, {
    "myworkflow.x.x": "101"
  }
]

--f3fd038395644de596c460257626edd8
Content-Disposition: form-data; name="workflowOptions"; filename="options.json"

{
  "jes_gcs_root": "gs://myworkflow-dev/workflows"
}

--f3fd038395644de596c460257626edd8--
```


## GET /api/workflows/:version/query

This endpoint allows for querying workflows based on the following criteria:

* `name`
* `id`
* `status`
* `label`
* `start` (start datetime with mandatory offset)
* `end` (end datetime with mandatory offset)
* `page` (page of results)
* `pagesize` (# of results per page)

Names, ids, and statuses can be given multiple times to include
workflows with any of the specified names, ids, or statuses. When
multiple names are specified, any workflow matching one of the names
will be returned. The same is true for multiple ids or statuses. When
more than one label is specified, only workflows associated to all of
the given labels will be returned. 

When a combination of criteria are specified, for example querying by 
names and statuses, the results must return workflows that match one of 
the specified names and one of the statuses. Using page and pagesize will
enable server side pagination.

Valid statuses are `Submitted`, `Running`, `Aborting`, `Aborted`, `Failed`, and `Succeeded`.  `start` and `end` should
be in [ISO8601 datetime](https://en.wikipedia.org/wiki/ISO_8601) format with *mandatory offset* and `start` cannot be after `end`.

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/query?start=2015-11-01T00%3A00%3A00-04%3A00&end=2015-11-04T00%3A00%3A00-04%3A00&status=Failed&status=Succeeded&page=1&pagesize=10"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/query?start=2015-11-01T00%3A00%3A00-04%3A00&end=2015-11-04T00%3A00%3A00-04%3A00&status=Failed&status=Succeeded&page=1&pagesize=10"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 133
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
  "results": [
    {
      "name": "w",
      "id": "fdfa8482-e870-4528-b639-73514b0469b2",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:52.000-05:00",
      "start": "2015-11-01T07:38:57.000-05:00"
    },
    {
      "name": "hello",
      "id": "e69895b1-42ed-40e1-b42d-888532c49a0f",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:30.000-05:00",
      "start": "2015-11-01T07:38:58.000-05:00"
    },
    {
      "name": "crasher",
      "id": "ed44cce4-d21b-4c42-b76d-9d145e4d3607",
      "status": "Failed",
      "end": "2015-11-01T07:45:44.000-05:00",
      "start": "2015-11-01T07:38:59.000-05:00"
    }
  ]
}
```

Labels have to be queried in key and value pairs separated by a colon, i.e. `label-key:label-value`. For example, if a batch of workflows was submitted with the following labels JSON:
```
{
  "label-key-1" : "label-value-1",
  "label-key-2" : "label-value-2"
}
```

A request to query for succeeded workflows with both labels would be:

cURL:
```
$ curl "http://localhost:8000/api/workflows/v1/query?status=Succeeded&label=label-key-1:label-value-1&label=label-key-2:label-value-2
```

HTTPie:
```
$ http "http://localhost:8000/api/workflows/v1/query?status=Succeeded&label=label-key-1:label-value-1&label=label-key-2:label-value-2
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 608
Content-Type: application/json; charset=UTF-8
Date: Tue, 9 May 2017 20:24:33 GMT
Server: spray-can/1.3.3

{
    "results": [
        {
            "end": "2017-05-09T16:07:30.515-04:00", 
            "id": "83fc23d5-48d1-456e-997a-087e55cd2e06", 
            "name": "wf_hello", 
            "start": "2017-05-09T16:01:51.940-04:00", 
            "status": "Succeeded"
        }, 
        {
            "end": "2017-05-09T16:07:13.174-04:00", 
            "id": "7620a5c6-a5c6-466c-994b-dd8dca917b9b", 
            "name": "wf_goodbye", 
            "start": "2017-05-09T16:01:51.939-04:00", 
            "status": "Succeeded"
        }
    ]
}
```

Query data is refreshed from raw data periodically according to the configuration value `services.MetadataService.metadata-summary-refresh-interval`.
This interval represents the duration between the end of one summary refresh sweep and the beginning of the next sweep.  If not specified the
refresh interval will default to 2 seconds.  To turn off metadata summary refresh, specify an infinite refresh interval value with "Inf".

```
services {
  MetadataService {
    config {
      metadata-summary-refresh-interval = "10 seconds"
    }
  }
}
```

## POST /api/workflows/:version/query

This endpoint allows for querying workflows based on the same criteria
as [GET /api/workflows/:version/query](#get-apiworkflowsversionquery).

Instead of specifying query parameters in the URL, the parameters
must be sent via the POST body. The request content type must be
`application/json`. The json should be a list of objects. Each json
object should contain a different criterion.

cURL:

```
$ curl -X POST --header "Content-Type: application/json" -d "[{\"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"page\": \"1\"}, {\"pagesize\": \"10\"}]" "http://localhost:8000/api/workflows/v1/query"
```

HTTPie:

```
$ echo "[{\"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"page\": \"1\"}, {\"pagesize\": \"10\"}]" | http "http://localhost:8000/api/workflows/v1/query"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 133
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
  "results": [
    {
      "name": "w",
      "id": "fdfa8482-e870-4528-b639-73514b0469b2",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:52.000-05:00",
      "start": "2015-11-01T07:38:57.000-05:00"
    },
    {
      "name": "hello",
      "id": "e69895b1-42ed-40e1-b42d-888532c49a0f",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:30.000-05:00",
      "start": "2015-11-01T07:38:58.000-05:00"
    },
    {
      "name": "crasher",
      "id": "ed44cce4-d21b-4c42-b76d-9d145e4d3607",
      "status": "Failed",
      "end": "2015-11-01T07:45:44.000-05:00",
      "start": "2015-11-01T07:38:59.000-05:00"
    }
  ]
}
```

## PATCH /api/workflows/:version/:id/labels

This endpoint is used to update multiple labels for an existing workflow. When supplying a label with a key unique to the workflow submission, a new label key/value entry is appended to that workflow's metadata. When supplying a label with a key that is already associated to the workflow submission, the original label value is updated with the new value for that workflow's metadata.

The [labels](#labels) must be a mapping of key/value pairs in JSON format that are sent via the PATCH body. The request content type must be
`application/json`.

cURL:

```
$ curl -X PATCH --header "Content-Type: application/json" -d "{\"label-key-1\":\"label-value-1\", \"label-key-2\": \"label-value-2\"}" "http://localhost:8000/api/workflows/v1/c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5/labels"
```

HTTPie:

```
$ echo '{"label-key-1":"label-value-1", "label-key-2": "label-value-2"}' | http PATCH "http://localhost:8000/api/workflows/v1/c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5/labels"
```

Response:
```
{ "id": "c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5",
  "labels":
    {
      "label-key-1": "label-value-1",
      "label-key-2": "label-value-2"
    }
}
```


## GET /api/workflows/:version/:id/status

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 74
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
    "id": "69d1d92f-3895-4a7b-880a-82535e9a096e",
    "status": "Succeeded"
}
```

## GET /api/workflows/:version/:id/outputs

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 241
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1",
    "outputs": {
        "three_step.cgrep.count": 8,
        "three_step.ps.procs": "/var/folders/kg/c7vgxnn902lc3qvc2z2g81s89xhzdz/T/stdout2814345504446060277.tmp",
        "three_step.wc.count": 8
    }
}
```

## GET /api/workflows/:version/:id/timing

This endpoint is meant to be used in a web browser.  It will show a Gantt Chart of a particular workflow.  The bars in the chart represent start and end times for individual task invocations.

![Timing diagram](http://i.imgur.com/EOE2HoL.png)

## GET /api/workflows/:version/:id/logs

This will return paths to the standard out and standard error files that were generated during the execution of all calls in a workflow.

A call has one or more standard out and standard error logs, depending on if the call was scattered or not. In the latter case, one log is provided for each instance of the call that has been run.

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 379
Content-Type: application/json; charset=UTF-8
Date: Mon, 03 Aug 2015 17:11:28 GMT
Server: spray-can/1.3.3

{
    "id": "b3e45584-9450-4e73-9523-fc3ccf749848",
    "logs": {
        "call.ps": [
            {
                "stderr": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-ps/stderr6126967977036995110.tmp",
                "stdout": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-ps/stdout6128485235785447571.tmp"
            }
        ],
        "call.cgrep": [
            {
                "stderr": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-cgrep/stderr6126967977036995110.tmp",
                "stdout": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-cgrep/stdout6128485235785447571.tmp"
            }
        ],
        "call.wc": [
            {
                "stderr": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-wc/stderr6126967977036995110.tmp",
                "stdout": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-wc/stdout6128485235785447571.tmp"
            }
        ]
    }
}
```

## GET /api/workflows/:version/:id/metadata

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

## POST /api/workflows/:version/:id/abort

cURL:

```
$ curl -X POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/abort
```

HTTPie:

```
$ http POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/abort
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 241
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1",
    "status": "Aborted"
}
```

## GET /api/workflows/:version/backends

This endpoint returns a list of the backends supported by the server as well as the default backend.

cURL:
```
$ curl http://localhost:8000/api/workflows/v1/backends
```

HTTPie:
```
$ http http://localhost:8000/api/workflows/v1/backends
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 379
Content-Type: application/json; charset=UTF-8
Date: Mon, 03 Aug 2015 17:11:28 GMT
Server: spray-can/1.3.3

{
  "supportedBackends": ["JES", "LSF", "Local", "SGE"],
  "defaultBackend": "Local"
}
```

## GET /api/workflows/:version/callcaching/diff

**Disclaimer**: This endpoint depends on hash values being published to the metadata, which only happens as of Cromwell 28.
Workflows run with prior versions of Cromwell cannot be used with this endpoint.
A `404 NotFound` will be returned when trying to use this endpoint if either workflow has been run on a prior version.

This endpoint returns the hash differences between 2 *completed* (successfully or not) calls.
The following query parameters are supported:

| Parameter |                                        Description                                        | Required |
|:---------:|:-----------------------------------------------------------------------------------------:|:--------:|
| workflowA | Workflow ID of the first call                                                             | yes      |
| callA     | Fully qualified name of the first call. **Including workflow name**. (see example below)  | yes      |
| indexA    | Shard index of the first call                                                             | depends  |
| workflowB | Workflow ID of the second call                                                            | yes      |
| callB     | Fully qualified name of the second call. **Including workflow name**. (see example below) | yes      |
| indexB    | Shard index of the second call                                                            | depends  |

About the `indexX` parameters: It is required if the call was in a scatter. Otherwise it should *not* be specified.
If an index parameter is wrongly specified, the call will not be found and the request will result in a 404 response.

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/callcaching/diff?workflowA=85174842-4a44-4355-a3a9-3a711ce556f1&callA=wf_hello.hello&workflowB=7479f8a8-efa4-46e4-af0d-802addc66e5d&callB=wf_hello.hello"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/callcaching/diff?workflowA=85174842-4a44-4355-a3a9-3a711ce556f1&callA=wf_hello.hello&workflowB=7479f8a8-efa4-46e4-af0d-802addc66e5d&callB=wf_hello.hello"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 1274
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 16:44:33 GMT
Server: spray-can/1.3.3

{
  "callA": {
    "executionStatus": "Done",
    "workflowId": "85174842-4a44-4355-a3a9-3a711ce556f1",
    "callFqn": "wf_hello.hello",
    "jobIndex": -1,
    "allowResultReuse": true
  },
  "callB": {
    "executionStatus": "Done",
    "workflowId": "7479f8a8-efa4-46e4-af0d-802addc66e5d",
    "callFqn": "wf_hello.hello",
    "jobIndex": -1,
    "allowResultReuse": true
  },
  "hashDifferential": [
    {
      "hashKey": "command template",
      "callA": "4EAADE3CD5D558C5A6CFA4FD101A1486",
      "callB": "3C7A0CA3D7A863A486DBF3F7005D4C95"
    },
    {
      "hashKey": "input count",
      "callA": "C4CA4238A0B923820DCC509A6F75849B",
      "callB": "C81E728D9D4C2F636F067F89CC14862C"
    },
    {
      "hashKey": "input: String addressee",
      "callA": "D4CC65CB9B5F22D8A762532CED87FE8D",
      "callB": "7235E005510D99CB4D5988B21AC97B6D"
    },
    {
      "hashKey": "input: String addressee2",
      "callA": "116C7E36B4AE3EAFD07FA4C536CE092F",
      "callB": null
    }
  ]
}
```

The response is a JSON object with 3 fields:

- `callA` reports information about the first call, including its `allowResultReuse` value that will be used to determine whether or not this call can be cached to.
- `callB` reports information about the second call, including its `allowResultReuse` value that will be used to determine whether or not this call can be cached to.
- `hashDifferential` is an array in which each element represents a difference between the hashes of `callA` and `callB`.

*If this array is empty, `callA` and `callB` have the same hashes*.

Differences can be of 3 kinds:

- `callA` and `callB` both have the same hash key but their values are different.
For instance, in the example above, 

```json
{
  "hashKey": "input: String addressee",
  "callA": "D4CC65CB9B5F22D8A762532CED87FE8D",
  "callB": "7235E005510D99CB4D5988B21AC97B6D"
}
```

indicates that both `callA` and `callB` have a `String` input called `addressee`, but different values were used at runtime, resulting in different MD5 hashes.

- `callA` has a hash key that `callB` doesn't have
For instance, in the example above, 

```json
{
  "hashKey": "input: String addressee2",
  "callA": "116C7E36B4AE3EAFD07FA4C536CE092F",
  "callB": null
}
```

indicates that `callA` has a `String` input called `addressee2` that doesn't exist in `callB`. For that reason the value of the second field is `null`.

- `callB` has a hash key that `callA` doesn't have. This is the same case as above but reversed.

If no cache entry for `callA` or `callB` can be found, the response will be in the following format:

```
HTTP/1.1 404 NotFound
Content-Length: 178
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 17:02:15 GMT
Server: spray-can/1.3.3

{
  "status": "error",
  "message": "Cannot find call 479f8a8-efa4-46e4-af0d-802addc66e5d:wf_hello.hello:-1"
}
```

If neither `callA` nor `callB` can be found, the response will be in the following format:


```
HTTP/1.1 404 NotFound
Content-Length: 178
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 17:02:15 GMT
Server: spray-can/1.3.3

{
  "status": "error",
  "message": "Cannot find calls 5174842-4a44-4355-a3a9-3a711ce556f1:wf_hello.hello:-1, 479f8a8-efa4-46e4-af0d-802addc66e5d:wf_hello.hello:-1"
}
```

If the query is malformed and required parameters are missing, the response will be in the following format:

```
HTTP/1.1 400 BadRequest
Content-Length: 178
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 17:02:15 GMT
Server: spray-can/1.3.3
{
  "status": "fail",
  "message": "Wrong parameters for call cache diff query:\nmissing workflowA query parameter\nmissing callB query parameter",
  "errors": [
    "missing workflowA query parameter",
    "missing callB query parameter"
  ]
}
```

## GET /engine/:version/stats

This endpoint returns some basic statistics on the current state of the engine. At the moment that includes the number of running workflows and the number of active jobs. 

cURL:
```
$ curl http://localhost:8000/engine/v1/stats
```

HTTPie:
```
$ http http://localhost:8000/engine/v1/stats
```

Response:
```
"date": "Sun, 18 Sep 2016 14:38:11 GMT",
"server": "spray-can/1.3.3",
"content-length": "33",
"content-type": "application/json; charset=UTF-8"

{
  "workflows": 3,
  "jobs": 10
}
```

## GET /engine/:version/version

This endpoint returns the version of the Cromwell engine.

cURL:
```
$ curl http://localhost:8000/engine/v1/version
```

HTTPie:
```
$ http http://localhost:8000/engine/v1/version
```

Response:
```
"date": "Sun, 18 Sep 2016 14:38:11 GMT",
"server": "spray-can/1.3.3",
"content-length": "33",
"content-type": "application/json; charset=UTF-8"

{
  "cromwell": 23-8be799a-SNAP
}
```

## GET /engine/:version/status

Cromwell will track the status of underlying systems on which it depends, typically connectivity to the database and to Docker Hub. The `status` endpoint will return the current status of these systems. 

Response:

This endpoint will return an `Internal Server Error` if any systems are marked as failing or `OK` otherwise. The exact response will vary based on what is being monitored but adheres to the following pattern. Each system has a boolean `ok` field and an optional array field named `messages` which will be populated with any known errors if the `ok` status is `false`:

```
{
  "System Name 1": {
    "ok": false,
    "messages": [
      "Unknown status"
    ]
  },
  "System Name 2": {
    "ok": true
  }
}

```
 


## Error handling
Requests that Cromwell can't process return a failure in the form of a JSON response respecting the following JSON schema:

```
{
  "$schema": "http://json-schema.org/draft-04/schema#",
  "description": "Error response schema",
  "type": "object",
  "properties": {
    "status": {
      "enum": [ "fail", "error"]
    },
    "message": {
      "type": "string"
    },
    "errors": {
      "type": "array",
      "minItems": 1,
      "items": { "type": "string" },
      "uniqueItems": true
    }
  },
  "required": ["status", "message"]
}
```

The `status` field can take two values:
> "fail" means that the request was invalid and/or data validation failed. "fail" status is most likely returned with a 4xx HTTP Status code.
e.g.

```
{
  "status": "fail",
  "message": "Workflow input processing failed.",
  "errors": [
    "Required workflow input 'helloworld.input' not specified."
  ]
}
```

> "error" means that an error occurred while processing the request. "error" status is most likely returned with a 5xx HTTP Status code.
e.g.

```
{
  "status": "error",
  "message": "Connection to the database failed."
}
```

The `message` field contains a short description of the error.

The `errors` field is optional and may contain additional information about why the request failed.