For built-in documentation of Cromwell command line usage, run the Cromwell JAR file with no arguments:

```
$ java -jar cromwell-<versionNumber>.jar
```

For example, `$ java -jar cromwell-29.jar`. You will get a usage message like the following:

```
cromwell 29
Usage: java -jar /path/to/cromwell.jar [server|run] [options] <args>...

  --help                   Cromwell - Workflow Execution Engine
  --version                
Command: server
Starts a web server on port 8000.  See the web server documentation for more details about the API endpoints.
Command: run [options] workflow-source
Run the workflow and print out the outputs in JSON format.
  workflow-source          Workflow source file.
  -i, --inputs <value>     Workflow inputs file.
  -o, --options <value>    Workflow options file.
  -t, --type <value>       Workflow type.
  -v, --type-version <value>
                           Workflow type version.
  -l, --labels <value>     Workflow labels file.
  -p, --imports <value>    A directory or zipfile to search for workflow imports.
  -m, --metadata-output <value>
                           An optional directory path to output metadata.
```

## --version

The `--version` option prints the version of Cromwell and exits.

## --help

The `--help` option prints the full help text above and exits.

## server

The `server` command runs Cromwell as a web server.  No arguments are accepted.
See the documentation for Cromwell's REST endpoints [here](#rest-api).

## run

The `run` command executes a single workflow in Cromwell.

### workflow-source
The `run` command requires a single argument for the workflow source file.
 
### --inputs
An optional file of workflow inputs.  Although optional, it is a best practice to use an inputs file to satisfy workflow
requirements rather than hardcoding inputs directly into a workflow source file.

### --options
An optional file of workflow options.  Some options are global (supported by all backends), while others are backend-specific.
See the [workflow options](#workflow-options) documentation for more details.

### --type
An optional parameter to specify the language for the workflow source.  Any value specified for this parameter is currently
ignored and internally the value `WDL` is used.

### --type-version
An optional parameter to specify the version of the language for the workflow source.  Currently any specified value is ignored.

### --labels
An optional parameter to specify a file of JSON key-value label pairs to associate with the workflow.

### --imports
You have the option of importing WDL workflows or tasks to use within your workflow, known as sub-workflows.
If you use sub-workflows within your primary workflow then you must include a zip file with the WDL import files.

For example, say you have a directory of WDL files:

```
wdl_library
└──cgrep.wdl
└──ps.wdl
└──wc.wdl
```

If you zip that directory into `wdl_library.zip`, then you can reference and use these WDLs within your primary WDL.

This could be your primary WDL:

```
import "ps.wdl" as ps
import "cgrep.wdl"
import "wc.wdl" as wordCount

workflow my_wf {

call ps.ps as getStatus
call cgrep.cgrep { input: str = getStatus.x }
call wordCount { input: str = ... }

}
```

Then to run this WDL without any inputs, workflow options, or metadata files, you would enter:

`$ java -jar cromwell-<versionNumber>.jar run my_wf.wdl --imports /path/to/wdl_library.zip`

### --metadata-output

You can include a path where Cromwell will write the workflow metadata JSON, such as start/end timestamps, status, inputs, and outputs. By default, Cromwell does not write workflow metadata.

This example includes a metadata path called `/path/to/my_wf.metadata`:

```
$ java -jar cromwell-<versionNumber>.jar run my_wf.wdl --metadata-output /path/to/my_wf.metadata
```

Again, Cromwell is very verbose. Here is the metadata output in my_wf.metadata:

```
{
  "workflowName": "my_wf",
  "submittedFiles": {
    "inputs": "{\"my_wf.hello.addressee\":\"m'Lord\"}",
    "workflow": "\ntask hello {\n  String addressee\n  command {\n    echo \"Hello ${addressee}!\"\n  }\n  output {\n    String salutation = read_string(stdout())\n  }\n  runtime {\n   
\n  }\n}\n\nworkflow my_wf {\n  call hello\n  output {\n     hello.salutation\n  }\n}\n",
    "options": "{\n\n}"
  },
  "calls": {
    "my_wf.hello": [
      {
        "executionStatus": "Done",
        "stdout": "/Users/jdoe/Documents/cromwell-executions/my_wf/cd0fe94a-984e-4a19-ab4c-8f7f07038068/call-hello/execution/stdout",
        "backendStatus": "Done",
        "shardIndex": -1,
        "outputs": {
          "salutation": "Hello m'Lord!"
        },
        "runtimeAttributes": {
          "continueOnReturnCode": "0",
          "failOnStderr": "false"
        },
        "callCaching": {
          "allowResultReuse": false,
          "effectiveCallCachingMode": "CallCachingOff"
        },
        "inputs": {
          "addressee": "m'Lord"
        },
        "returnCode": 0,
        "jobId": "28955",
        "backend": "Local",
        "end": "2017-04-19T10:53:25.045-04:00",
        "stderr": "/Users/jdoe/Documents/cromwell-executions/my_wf/cd0fe94a-984e-4a19-ab4c-8f7f07038068/call-hello/execution/stderr",
        "callRoot": "/Users/jdoe/Documents/cromwell-executions/my_wf/cd0fe94a-984e-4a19-ab4c-8f7f07038068/call-hello",
        "attempt": 1,
        "executionEvents": [
          {
            "startTime": "2017-04-19T10:53:23.570-04:00",
            "description": "PreparingJob",
            "endTime": "2017-04-19T10:53:23.573-04:00"
          },
          {
            "startTime": "2017-04-19T10:53:23.569-04:00",
            "description": "Pending",
            "endTime": "2017-04-19T10:53:23.570-04:00"
          },
          {
            "startTime": "2017-04-19T10:53:25.040-04:00",
            "description": "UpdatingJobStore",
            "endTime": "2017-04-19T10:53:25.045-04:00"
          },
          {
            "startTime": "2017-04-19T10:53:23.570-04:00",
            "description": "RequestingExecutionToken",
            "endTime": "2017-04-19T10:53:23.570-04:00"
          },
          {
            "startTime": "2017-04-19T10:53:23.573-04:00",
            "description": "RunningJob",
            "endTime": "2017-04-19T10:53:25.040-04:00"
          }
        ],
        "start": "2017-04-19T10:53:23.569-04:00"
      }
    ]
  },
  "outputs": {
    "my_wf.hello.salutation": "Hello m'Lord!"
  },
  "workflowRoot": "/Users/jdoe/Documents/cromwell-executions/my_wf/cd0fe94a-984e-4a19-ab4c-8f7f07038068",
  "id": "cd0fe94a-984e-4a19-ab4c-8f7f07038068",
  "inputs": {
    "my_wf.hello.addressee": "m'Lord"
  },
  "submission": "2017-04-19T10:53:19.565-04:00",
  "status": "Succeeded",
  "end": "2017-04-19T10:53:25.063-04:00",
  "start": "2017-04-19T10:53:23.535-04:00"
}
```