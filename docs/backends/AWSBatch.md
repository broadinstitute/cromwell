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

To allow Cromwell to talk to AWS, we have a few authentication schemes that we can specify in the configuration that Cromwell can accept:
-  `default`: This uses the [default authentication provider](https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html), which uses the standard AWS search paths:
    - Environment Variables - `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`
    - Java system properties - `aws.accessKeyId` and `aws.secretKey`
    - Default credential profiles file - Created by the AWS CLI, typically located at `~/.aws/credentials`
    - _Instance profile credentials_ - Only relevant on EC2 instances
- `custom_keys`: This requires the definition of "access-key" and "secret-key"
- `assume_role`: This requires the definition of "base-auth" and "role-arn". An optional "external-id" can be provided.

Specifically, within the `aws.auths` array, you can provide a list of authentication options:

```hocon

# ... within the aws.conf
aws {
    application-name = "cromwell"
    auths = [
        { name="default", scheme="default" },
        { name="custom_keys", access-key="<my-access-key>", role-arn="<my-role-arn>" },
        { name="assume_role", base-auth="<my-base-auth>", role-arn="<my-role-arn>"}
    ]
}
```

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
