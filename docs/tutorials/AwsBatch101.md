## Getting started on AWS with AWS Batch (beta)

### Prerequisites

This tutorial page relies on completing the previous tutorial:

* [Downloading Prerequisites](FiveMinuteIntro.md)

### Goals

At the end of this tutorial you'll have configured your local environment to run workflows using Cromwell on AWS Batch.

### Let's get started!

To create all the resources for running a Cromwell server on AWS using CloudFormation, launch the [Cromwell Full Stack Deployment](https://docs.opendata.aws/genomics-workflows/orchestration/cromwell/cromwell-overview/).  Alternatively, this page will walk through the specific steps to configure and run a local Cromwell server using AWS Batch.

1. [Authenticating a local Cromwell server with AWS](#authenticating-a-local-cromwell-server-with-aws)
2. [Configuring the AWS environment](#configuring-the-aws-environment)
3. [Configuring Cromwell](#configuring-cromwell)
4. [Workflow Source Files](#workflow-source-files)
5. [Running Cromwell and AWS](#running-cromwell-and-aws)
6. [Outputs](#outputs)

#### Authenticating a local Cromwell server with AWS

The easiest way to allow a local Cromwell server to talk to AWS is to:

1. Install the AWS CLI through Amazon's [user guide](https://docs.aws.amazon.com/cli/latest/userguide/installing.html).
2. [Configure the AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html) by calling `aws configure` (provide your `Access Key` and `Secret Access Key` when prompted).

Cromwell can access these credentials through the default authentication provider. For more options, see the [Configuring authentication of Cromwell with AWS](/backends/AWSBatch#configuring-authentication) section below.


#### Configuring the AWS environment

Next you'll need the following setup in your AWS account:
- The core set of resources (S3 Bucket, IAM Roles, AWS Batch)
- Custom Compute Resource (Launch Template or AMI) with Cromwell Additions

Information and instructions to setup an AWS environment to work properly with Cromwell can be found on [AWS for Genomics Workflow](https://docs.opendata.aws/genomics-workflows/core-env/introduction/). By deploying the CloudFormation templates provided by AWS, the stack will output the S3 bucket name and two AWS Batch queue ARNs (default and high-priority) used in the Cromwell configuration.




#### Configuring Cromwell

Now we're going to configure Cromwell to use the AWS resources we just created by updating a `*.conf` file to use the `AWSBackend` at runtime. This requires three pieces of information:

- The [AWS Region](https://docs.aws.amazon.com/general/latest/gr/rande.html) where your resources are deployed.
- S3 bucket name where Cromwell will store its execution files.
- The ARN of the AWS Batch queue you want to use for your tasks.

You can replace the placeholders (`<your region>`, `<your-s3-bucket-name>` and `<your-queue-arn>`) in the following config:

##### `aws.conf`

```hocon
include required(classpath("application"))

aws {

  application-name = "cromwell"
  auths = [
    {
      name = "default"
      scheme = "default"
    }
  ]
  region = "<your-region>"
}

engine {
  filesystems {
    s3.auth = "default"
  }
}

backend {
  default = "AWSBatch"
  providers {
    AWSBatch {
      actor-factory = "cromwell.backend.impl.aws.AwsBatchBackendLifecycleActorFactory"
      config {
        
        numSubmitAttempts = 6
        numCreateDefinitionAttempts = 6

        // Base bucket for workflow executions
        root = "s3://<your-s3-bucket-name>/cromwell-execution"

        // A reference to an auth defined in the `aws` stanza at the top.  This auth is used to create
        // Jobs and manipulate auth JSONs.
        auth = "default"

        default-runtime-attributes {
          queueArn: "<your arn here>"
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

For more information about this configuration or how to change the behaviour of AWS Batch, visit the [AWS Backend](/backends/AWSBatch) page.


#### Workflow Source Files 

Lastly, create an example workflow to run. We're going to define a simple workflow that will `echo` a string to the console and return the result to Cromwell. Within AWS Batch (like other cloud providers), we're required to specify a Docker container for every task.

##### `hello.wdl`

```wdl
task hello {
  String addressee = "Cromwell"
  command {
    echo "Hello ${addressee}! Welcome to Cromwell . . . on AWS!"
  }
  output {
    String message = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  call hello
  output { hello.message }
}
```

#### Running Cromwell and AWS 

Provided all of the files are within the same directory, we can run our workflow with the following command:

> **Note**: You might have a different Cromwell version number here

```bash
java -Dconfig.file=aws.conf -jar cromwell-36.jar run hello.wdl
```

This will:
1. Start Cromwell in `run` mode,
2. Prepare `hello.wdl` as a job and submit this to your AWS Batch queue. You can monitor the job within your [AWS Batch dashboard](https://console.aws.amazon.com/batch/home).
3. Run the job, write execution files back to S3, and report progress back to Cromwell.

#### Outputs

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

You might find the following tutorials and guides interesting to tackle next:

* [Server Mode](/tutorials/ServerMode)
* [AWS Batch Backend](/backends/AWSBatch)
* [Persisting Data Between Restarts](/tutorials/PersistentServer)
