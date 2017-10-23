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