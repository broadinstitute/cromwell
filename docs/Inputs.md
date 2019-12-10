## Inputs

Inputs are specified as a JSON file. The JSON file can be passed to cromwell by using the
`-i` of `--inputs` command line flag. See the [command line documentation](CommandLine.md) 
for more information. Check the [API documentation](api/RESTAPI.md) for information on how
to pass inputs to a cromwell server.

### Example

Given the following workflow

***hello.wdl***
```WDL
version 1.0
task hello {
  input {
    String addressee
    String? optionalText
  }  
  command {
    echo "Hello ~{addressee}! Welcome to Cromwell . . . on Google Cloud!"
    ~{"echo " + optionalText}
  }
  output {
    String message = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  input {
      String addressee
  }
  call hello {
    input:
      addressee=addressee
  }

  output {
     String message = hello.message
  }
}
```

The following inputs file contains all the required inputs.
***hello.inputs***
```JSON
{
  "wf_hello.addressee": "World"
}
```

The inputs file can also be used to specify task level inputs.

```json
{
  "wf_hello.addressee": "World",
  "wf_hello.hello.optionalText": "This is submitted into the task input directly."
}
```

### Advanced: namespacing in inputs files.

Cromwell uses namespacing to make sure all inputs are addressable in a comprehendable manner.
In the above example the `wf_hello` workflow has its own namespace. In this namespace reside:
+ `addressee`, an input string
+ `hello`, a task
+ `message`, an output string

In the inputs file all inputs are addressable by namespace. So `addressee` is accessible as
`wf_hello.addressee`. This also allows access to all inputs of the task `hello` such as
`wf_hello.hello.optionalText`.

If a more complicated workflow with lots of subworkflows is used the namespaces can be used as
well.

```wdl
version 1.0 

import "subworkflow1.wdl" as wf_one
import "subworkflow2.wdl" as wf_two

workflow MainWorkflow {
  input {
    File dataFile
    String method      
  }    
  call wf_one.Workflow as WorkflowOne {
    input:
       data = dataFile
  }                    
  call wf_two.Workflow as WorkflowTwo {
     input:            
       data = WorkflowOne.analysedData,
       method = method
  }        
  output {
     File processedData = WorkflowTwo.methodicallyAnalyzedData
     File report = WorkflowTwo.report            
  }
}
```

The following inputs are required
```json
{
  "MainWorkflow.dataFile": "/path/to/data",
  "MainWorkflow.method": "thorough"      
}
```

Using namespaces, subworkflow tasks can also be addressed.
```json
{
  "MainWorkflow.dataFile": "/path/to/data",
  "MainWorkflow.method": "thorough",
  "MainWorkflow.WorkflowOne.nested_task.obscureParameter": 42,
  "MainWorkflow.WorkflowTwo.summarize.reportType": "html"
}
```

This feature allows users to specify inputs for deeply nested subworkflows as well:
```json
{
  "MainWorkflow.dataFile": "/path/to/data",
  "MainWorkflow.method": "thorough",
  "MainWorkflow.WorkflowOne.AnotherSubworkflow.DeeplyNestedWorkflow.analysisTask.kmerSize": 22
}
```