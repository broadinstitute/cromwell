## Getting started on Alibaba Cloud with the Batch Compute Service

### Prerequisites

This tutorial page relies on completing the previous tutorials:

- [Configuration Files](ConfigurationFiles.md)

### Goals

In this tutorial you'll learn to run the first workflow against the Batch Compute service on Alibaba Cloud.

### Let's get started!

####

#### Configuring Alibaba Cloud

- Go to <a href="https://www.aliyun.com/" target="_blank">Alibaba Cloud</a> and activate <a href="https://www.aliyun.com/product/oss">Alibaba Cloud OSS</a> and <a href="https://www.aliyun.com/product/batchcompute">Alibaba Cloud BatchCompute</a> services. 
- Follow <a href="https://help.aliyun.com/document_detail/63724.html" target="_blank">AccessKey Guide</a> to retrieve an access-id and access-key pair. We will refer to this pair as `<test-access-id>` and `<test-access-key>`, respectively.
- Log on to the <a href="https://oss.console.aliyun.com/" target="_blank">OSS console</a> and choose a region to create a new bucket. We will use `<test-region>` and `<test-bucket>` to refer the chosen region and bucket.
- Find the corresponding OSS endpoint in <a href="https://help.aliyun.com/document_detail/31837.html" target="_blank">OSS region and endpoint</a>. We will refer to it as `<test-oss-endpoint>`.

#### Preparing workflow source files

Copy over the sample `echo.wdl` and `echo.inputs` files to the same directory as the Cromwell jar. 
This workflow takes a string value as an output file name and writes "Hello World!" to the file. 

***echo.wdl***

```
task echo {
  String out

  command {
    echo Hello World! > ${out}
  }

  output {
    File outFile = "${out}"
    Array[String] content = read_lines(outFile)
  }
}

workflow wf_echo {
  call echo
  output {
    echo.outFile
    echo.content
  }
}
```

***echo.inputs***

```
{
  "wf_echo.echo.out": "output"
}
```

#### Configuration file for Alibaba Cloud

Copy over the sample `bcs.conf` file to the same directory that contains your sample WDL, inputs and the Cromwell jar. Replace `<test-bucket>`, `<test-region>`, `<test-access-id>`, `<test-access-key>`, `<test-oss-endpoint>` in the configuration file with actual values.  

***bcs.conf***

```
include required(classpath("application"))

backend {
  default = "BCS"
  
  providers {
    BCS {
      actor-factory = "cromwell.backend.impl.bcs.BcsBackendLifecycleActorFactory"
      config {
        root = "oss://<test-bucket>/cromwell-dir"
        region = "<test-region>"
        access-id = "<test-access-id>"
        access-key = "<test-access-key>"
        
        filesystems {
          oss {
            auth {
              endpoint = "<test-oss-endpoint>"
              access-id = "<test-access-id>"
              access-key = "<test-access-key>"
            }
          }
        }
        
        default-runtime-attributes {
          failOnStderr: false
          continueOnReturnCode: 0
          cluster: "OnDemand ecs.sn1ne.large img-ubuntu"
          vpc: "192.168.0.0/16"
        } 
      }
    }
  }
}
```

#### Run workflow

`java -Dconfig.file=bcs.conf -jar cromwell.jar run echo.wdl --inputs echo.inputs`

#### Outputs

The end of your workflow logs should report the workflow outputs. 

```
[info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "wf_echo.echo.outFile": "oss://<test-bucket>/cromwell-dir/wf_echo/38b088b2-5131-4ea0-a161-4cf2ca8d15ac/call-echo/output",
    "wf_echo.echo.content": ["Hello World!"]
  },
  "id": "38b088b2-5131-4ea0-a161-4cf2ca8d15ac"
}
```
