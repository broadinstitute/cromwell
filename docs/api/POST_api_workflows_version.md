_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/versions page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---



This endpoint accepts a POST request with a `multipart/form-data` encoded body.  The form fields that may be included are:

* `workflowSource` - *Required* Contains the workflow source file to submit for execution.
* `workflowType` - *Optional* The type of the `workflowSource`.
   If not specified, returns the `workflow-options.default.workflow-type` configuration value with a default of `WDL`.
* `workflowTypeVersion` - *Optional* The version of the `workflowType`, for example "draft-2".
   If not specified, returns the `workflow-options.default.workflow-type-version` configuration value with no default.
* `workflowInputs` - *Optional* JSON file containing the inputs.  For WDL workflows a skeleton file can be generated from [wdltool](https://github.com/broadinstitute/wdltool) using the "inputs" subcommand.
* `workflowInputs_n` - *Optional* Where `n` is an integer. JSON file containing the 'n'th set of auxiliary inputs.
* `workflowOptions` - *Optional* JSON file containing options for this workflow execution.  See the [run](#run) CLI sub-command for some more information about this.
* `customLabels` - *Optional* JSON file containing a set of custom labels to apply to this workflow. See [Labels](/labels) for the expected format.
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