## Getting started with AWS Batch

### Prerequisites

This tutorial page relies on completing the previous tutorial:

* [Downloading Prerequisites](FiveMinuteIntro.md)

### Goals

At the end of this tutorial you'll have run your first workflow against AWS Batch

### Let's get started!


**Configuring the AWS environment**

Install the <a href="https://docs.aws.amazon.com/cli/latest/userguide/installing.html" target="_blank">AWS CLI</a>.

<a href="https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html" target="_blank">Configure the CLI for use.</a>
Note that if using an EC2 instance, the <a href="https://docs.aws.amazon.com/cli/latest/userguide/cli-metadata.html" target="_blank">
instance metadata, provided by an instance role</a> is generally preferred.
Also, any EC2 instance using <a href="https://aws.amazon.com/amazon-linux-ami/" target="_blank">Amazon Linux</a>
comes with the AWS CLI pre-installed.

Create a <a href="https://docs.aws.amazon.com/AmazonS3/latest/dev/UsingBucket.html">S3 bucket</a> to hold Cromwell execution directories.
We will refer to this bucket as `s3-bucket-name`, and the full identifier as `s3://s3-bucket-name`.
`aws s3 mb s3://<s3-bucket-name>` 


**Workflow Source Files**

Copy over the sample `hello.wdl` and `hello.inputs` files to the same directory as the Cromwell jar. 
This workflow takes a string value as specified in the inputs file and writes it to stdout. 


***hello.wdl***
```
task hello {
  String addressee
  command {
    echo "Hello ${addressee}! Welcome to Cromwell . . . on AWS!"
  }
  output {
    String message = read_string(stdout())
  }
  runtime {
    docker: "alpine:latest"
  }
}

workflow wf_hello {
  call hello

  output {
     hello.message
  }
}
```

***hello.inputs***
```
{
  "wf_hello.hello.addressee": "World"
}
```

**AWS Configuration File**

Copy over the sample `aws.conf` file utilizing
<a href="https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html" target="_blank">the default credential provider chain</a>
to the same directory that contains your sample WDL, inputs and Cromwell jar.

Replace `<s3-bucket-name>`in the configuration file with the appropriate bucket name.

***aws.conf***
```
include required(classpath("application"))

aws {

  application-name = "cromwell"

  auths = [
    {
      name = "default"
      scheme = "default"
    }
  ]

  region = "default" // uses region from ~/.aws/config set by aws configure command,
                     // or us-east-1 by default
}

engine {
  filesystems {
    s3 {
      auth = "default"
    }
  }
}

backend {
  default = "AWSBATCH"
  providers {
    AWSBATCH {
      actor-factory = "cromwell.backend.impl.aws.BatchBackendLifecycleActorFactory"
      config {
        compute-service-account = "default"

        // Base bucket for workflow executions
        root = "s3://<s3-bucket-name>/cromwell-execution"

        // Polling for completion backs-off gradually for slower-running jobs.
        // This is the maximum polling interval (in seconds):
        maximum-polling-interval = 600

        // Optional Dockerhub Credentials. Can be used to access private docker images.
        dockerhub {
          // account = ""
          // token = ""
        }

        genomics {
          // A reference to an auth defined in the `aws` stanza at the top.  This auth is used to create
          // Pipelines and manipulate auth JSONs.
          auth = "application-default"
          // Endpoint for APIs, no reason to change this under normal circumstances.
          endpoint-url = ""
        }

        filesystems {
          s3 {
            // A reference to a potentially different auth for manipulating files via engine functions.
            auth = "default"
          }
        }
      }
    }
  }
}
```

**Run Workflow**

`java -Dconfig.file=aws.conf -jar cromwell-29.jar run hello.wdl -i hello.inputs`

**Outputs**

The end of your workflow logs should report the workflow outputs.

```
[info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "wf_hello.hello.message": "Hello World! Welcome to Cromwell . . . on AWS!"
  },
  "id": "08213b40-bcf5-470d-b8b7-1d1a9dccb10e"
}
```

Success!

### Next steps

You might find the following tutorials interesting to tackle next:

* [Persisting Data Between Restarts](PersistentServer)
* [Server Mode](ServerMode.md)
