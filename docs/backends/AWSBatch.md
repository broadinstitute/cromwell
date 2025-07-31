# AWS Batch Backend

AWS Batch is a set of batch management capabilities that dynamically provision the optimal quantity and type of compute resources (e.g., CPU or memory optimized instances) based on the volume and specific resource requirements of the batch jobs submitted.

This section provides details on how to configure the AWS Batch backend with Cromwell. For instructions on common configuration and deployment tutorial, see [Getting started with AWS Batch](https://cromwell.readthedocs.io/en/develop/tutorials/AwsBatch101/). 



## Resources and Runtime Attributes

Cromwell and AWS Batch recognizes number of runtime attributes, more information can be found in the [customize tasks](/RuntimeAttributes#recognized-runtime-attributes-and-backends) page.


## Running Cromwell on an EC2 instance

Cromwell can be run on an EC2 instance and submit jobs to AWS Batch, AWS provide [CloudFormation stacks and guides](https://docs.opendata.aws/genomics-workflows/orchestration/cromwell/cromwell-overview/) to building the correct IAM permissions.



## Scaling Requirements
For a Cromwell server that will run multiple workflows, or workflows with many steps (e.g. ones with large scatter steps), it is recommended to setup a database to store workflow metadata.  The application config file will expect a SQL database location. Follow [these instructions](https://docs.aws.amazon.com/AmazonRDS/latest/AuroraUserGuide/aurora-serverless.create.html) on how to create a serverless Amazon Aurora database. 

## Configuring Cromwell for AWS Batch

Within the `*.conf` file, you have a number of options to change the Cromwell's interaction with AWS Batch.

### Filesystems
> More information about filesystems can be found on the [Filesystems page](/filesystems/Filesystems/).
> 

Amazon's S3 storage is a supported filesystem in both the engine and backend, this means that S3 files can be referenced at a workflow level, and as input files, provided they are prefixed by `'s3://'`.

* filesystems
* filesystems.s3.auth
* filesystems.s3.caching.duplication-strategy

### Configuring Authentication

To allow Cromwell to talk to AWS, the `default` authentication scheme uses the [default authentication provider](https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html) with the following AWS search paths:
- Environment Variables - `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`
- Java system properties - `aws.accessKeyId` and `aws.secretKey`
- Default credential profiles file - Created by the AWS CLI, typically located at `~/.aws/credentials`
- _Instance profile credentials_ - Only relevant on EC2 instances

### Allowing private Docker containers

AWS Batch allows the use of private Docker containers by providing `dockerhub` credentials. Under the specific backend's configuration, you can provide the following object:

```hocon
(backend.providers.AWSBatch.config.)dockerhub = {
  // account = ""
  // token = ""
}
```

### More configuration options

* `(backend.providers.AWSBatch.config.)concurrent-job-limit` specifies the number of jobs that Cromwell will allow to be running in AWS at the same time. Tune this parameter based on how many nodes are in the compute environment.
* `(backend.providers.AWSBatch.config.)root` points to the S3 bucket where workflow outputs are stored. This becomes a path on the root instance, and by default is cromwell_root. This is monitored by preinstalled daemon that expands drive space on the host, ie AWS EBS autoscale.  This path is used as the 'local-disk' for containers.

## Workflow Options

The AWS Batch backend supports the following workflow options:

### aws_batch_job_role_arn

The `aws_batch_job_role_arn` workflow option allows you to specify an IAM role ARN that will be associated with the AWS Batch job's container. This enables the container to assume the specified role and access AWS resources according to the permissions granted to that role.

**Usage:**

Include this option in your workflow options JSON file:

```json
{
  "aws_batch_job_role_arn": "arn:aws:iam::123456789012:role/MyJobRole"
}
```

Or specify it on the command line:

```bash
java -jar cromwell.jar run workflow.wdl -o workflow_options.json
```

**Example use cases:**
- Accessing S3 buckets with specific permissions
- Invoking other AWS services (e.g., Lambda, SNS, SQS)
- Cross-account resource access

**Note:** The IAM role must have a trust relationship that allows the AWS Batch service to assume it. The role should include the following trust policy:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "ecs-tasks.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
```

### aws_batch_script_bucket_prefix

The `aws_batch_script_bucket_prefix` workflow option allows you to specify a custom prefix (subfolder) within the script bucket where execution scripts will be stored. This enables better organization and isolation of scripts for different workflows, projects, or teams.

**Usage:**

Include this option in your workflow options JSON file:

```json
{
  "aws_batch_script_bucket_prefix": "project-x/workflow-123/scripts"
}
```

**How it works:**
- By default (when not specified), scripts are stored in the `scripts/` folder within the configured `scriptBucketName`
- When you specify a prefix, scripts will be stored directly at that prefix location
- The prefix should include any subfolder structure you want (including `/scripts` if desired)
- If no prefix is specified or an empty string is provided, the default `scripts/` location is used
- A trailing slash will be added automatically if not present to ensure proper S3 key formation (e.g., `my-prefix` becomes `my-prefix/`)

**Example:**
If your `scriptBucketName` is configured as `my-cromwell-scripts` and you set:
```json
{
  "aws_batch_script_bucket_prefix": "team-alpha/project-genomics/scripts"
}
```

Scripts will be stored at: `s3://my-cromwell-scripts/team-alpha/project-genomics/scripts/[script-hash]`

Or for a custom structure without the `scripts` subfolder:
```json
{
  "aws_batch_script_bucket_prefix": "workflows/2024/genomics-pipeline"
}
```

Scripts will be stored at: `s3://my-cromwell-scripts/workflows/2024/genomics-pipeline/[script-hash]`

**Benefits:**
- **Organization**: Group scripts by workflow, project, team, or any other logical structure
- **Access Control**: Easier to implement S3 bucket policies based on prefixes
- **Cleanup**: Simpler to identify and remove old workflow scripts
- **Multi-tenancy**: Different teams or projects can have isolated script locations
