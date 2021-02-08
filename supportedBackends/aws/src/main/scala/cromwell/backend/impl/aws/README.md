AWS Batch Backend Architecture
==============================

Overview
--------

The architecture of the code base follows very closely to the Google version.
Probably a little too closely, and lots of code was lifted from the Google
backend originally, then modified to work with AWS.

Fundamentally, Google Pipelines API (a.k.a. PAPI) works pretty differently from
AWS Batch. In Pipelines, all the infrastructure is completely managed by Google,
while AWS Batch exposes that infrastructure to a large degree so that customers
can fine tune it as necessary. An implementation that uses Fargate might be an 
alternative that is closer or an implementation that uses Step Functions although
that would be a separate backend.

From a Cromwell perspective, this means that unlike Pipelines, where
infrastructure details are defined in the configuration or the WDL, in AWS
Batch, these configuration details are handled outside. All the AWS Batch
backend needs to know is "what is the ARN for the job Queue"?

A good example of the difference can be seen in the 'disks' configuration. In
Pipelines, you need to specify the type of disk and size. In AWS, this will
be defined instead when you setup your environment (more on that later), so
all the AWS backend really needs to know is what mount points you need
defined.

This infrastructure and all the associated configuration still exists; however,
it is moved out of the Cromwell configuration.

AWS Batch
---------

Because AWS Batch is so different from PAPI, those familiar only with PAPI
would be best off with an overview of AWS Batch. If you are familiar with
the workings of Batch, feel free to skip this section, and move on.

[AWS Batch](https://aws.amazon.com/batch/) fundamentally is a service to allow batch jobs to run easily and
efficiently. To use it effectively, however, you need to understand its own
technical stack. To create a job, you need a "Job Queue". That job queue allows
jobs to be scheduled onto one or more "Compute Environments". This can
be managed through AWS Batch, but when AWS Batch sets up a compute environment,
it's simply setting up an Elastic Container Service (ECS) Cluster. The ECS
cluster, in turn is just a few managed CloudFormation templates, that is
controlling an AutoScaling group of EC2 instances.

What really makes an ECS instance an ECS instance is the presence of a configured
[Amazon ECS agent](https://github.com/aws/amazon-ecs-agent). This agent polls
the ECS service to determine if there are any tasks to run. An AWS Batch Job
will be turned into an ECS task. From there, the agent will pick it up and
manage it, sending updates back to ECS (and from there AWS Batch) on a regular
basis.

There are some limits that will impact the design of the AWS Batch backend, and
will be discussed later. These are:

* [AWS Batch Limits](https://docs.aws.amazon.com/batch/latest/userguide/service_limits.html)
* [8k container overrides limit.](https://docs.aws.amazon.com/cli/latest/reference/ecs/run-task.html)

The ECS workers used by the AWS Batch backend can be any instance type and should
be based on an AMI running the ECS agent and docker. An ECS optimized AMI is recommended.
An EC2 LaunchTemplate is used to provide some additional "on first boot" configuration that:
1. Installs AWS CLI v2,
1. Installs a script to mount an EBS as a `btrfs` file system that will auto-expand,
1. Configures docker to use that file system so that the "filesystem" of the container
will auto-expand,
1. Installs a `fetch_and_run.sh` script that allows the container to download 
generated shell scripts from S3 that contain the instructions of the workflow
task 

```text
                  +-------------+
                  |             |
                  |  AWS Batch  |
                  |             |
                  +------+------+
                         |
                         |
                         |
                         |
                         |
        +----------------v------------------+
        |                                   |
        |  Elastic Container Service (ECS)  |
        |                                   |
        +----------------+------------------+
                         |
                         |
                         |
                         |
                         |
+------------------------v-------------------------+
|                                                  |
|  AutoScaling Group                               |
|                                                  |
| +---------------------------------+              |
| |                                 |              |
| |  EC2 Instance                   |              |
| |                                 |              |
| |  +--------------------+         |              |
| |  |                    |         |              |
| |  |  Docker Container  |         |              |
| |  |                    |         |              |
| |  +--------------------+  ...    |              |
| |                                 |              |
| +---------------------------------+     ...      |
|                                                  |
+--------------------------------------------------+

```

Cromwell AWS Batch Backend
--------------------------

There are several scala classes as part of the AWS Batch Backend, but
the primary classes involved in running the backend are shown below. The
arrows represent the flow of job submission.

```text
    +----------------------------------------+
    |                                        |
    |  AwsBatchBackendLifecycleActorFactory  |
    |                                        |
    +------------------+---------------------+
                       |
                       |
                       |
                       |
                       |
    +------------------v----------------------+
    |                                         |
    |  AwsBatchAsyncBackendJobExecutionActor  |
    |                                         |
    +------------------+----------------------+
                       |
                       |
                       |
                       |
                       |
               +-------v-------+                 +-------------------------+
               |               |                 |                         |
               |  AwsBatchJob  +----------------->  AwsBatchJobDefinition  |
               |               |                 |                         |
               +---------------+                 +-------------------------+
```

1. The `AwsBatchBackendLifecycleActorFactory` class is configured by the user
   as the Cromwell backend. This factory provides an object from the
   `AwsBatchAsyncBackendJobExecutionActor` class to create and manage the job.
2. The `AwsBatchAsyncBackendJobExecutionActor` creates and manages the job.
   The job itself is encapsulated by the functionality in `AwsBatchJob`.
3. `AwsBatchJob` is the primary interface to AWS Batch. It creates the
   necessary `AwsBatchJobDefinition`, then submits the job using the SubmitJob
   API.
4. `AwsBatchJobDefinition` is responsible for the creation of the job definition.
   In AWS Batch, every job must have a definition. Note that the job definition
   can be overridden by the `SubmitJob`, so the `JobDefinition` contains core information such
   as the docker image type while the `SubmitJob` contains details that are more related to
   the actual task.

AWS Batch Job Instantiation
---------------------------
```text
             +--------------------+
             |                    |
             |  Cromwell Backend  |
             |                    |
             +---------+----------+
                       |
                       |
                   SubmitJob
                       |
                       |
                +------v------+
                |             |
                |  AWS Batch  |
                |             |
                +------^------+
                       |
                       |
                     Polls
                       |
                       |
                +------+------+
                |             |
                |  ECS Agent  |
                |             |
                +------+------+
                       |
           Creates, Launches and Monitors
                       |
              +--------v---------+ 
              |                  |
              |  Task Container  |
              |                  |
              +------------------+

```

When a Cromwell task begins, the Cromwell backend will call the SubmitJob
API of AWS Batch. From there, the backend will call the AWS Batch `DescribeJobs`
API to provide status to the Cromwell engine as requested.

Once the job is Submitted in AWS Batch, one of the EC2 instances assigned
to the compute environment (a.k.a. ECS Cluster) with a running agent will
pick up the Cromwell Job/AWS Batch Job/ECS Task and run it. Importantly,
AWS Batch calls ECS' `RunTask` API when submitting the job. It uses the
task definition, and overrides both the command text and the environment
variables.

Input files are read into the container from S3 and output files are copied back to
S3. Three additional files are also written to the S3 bucket using the names of these
environment variables:

* AWS_CROMWELL_RC_FILE (the return code of the task)
* AWS_CROMWELL_STDOUT_FILE (STDOUT of the task)
* AWS_CROMWELL_STDERR_FILE (STDERR of the task)

These files are placed in the correct location in S3 after task execution. In addition
STDOUT and STDERR are fed to the tasks cloudwatch log.

Input and Command Compression
-----------------------------

NOTE: All limits in this section are subject to change

In testing, specifically with large fan-in operations such as the MergeVCFs
task of the Haplotype caller test, that the container overrides length limit
of 8k was being exceeded. There are several limits described on AWS Batch,
and a limit for container overrides on ECS, all of which should be considered.

* Maximum payload size for RegisterJobDefinition calls: 24KiB
* Maximum payload size for SubmitJob calls: 30KiB
* Maximum JSON payload for ECS RunTask containerOverrides values: 8KiB

Effective limits, however, are much, much smaller. While both AWS Batch and
ECS have command and environment as part of their Job/Task definitions
respectively, AWS Batch passes both command and environment through to ECS
based solely on RunTask. While a lot of effort was initially placed on
balancing payloads between RegisterJobDefinition and SubmitJob, because
everything is passed as an override to RunTask, we're gated by the ECS
RunTask 8KiB limit.

Dependencies
------------

Two dependencies were added to Cromwell as part of this project, though
existing Cromwell dependencies were also leveraged. These two are:

* AWS Java SDK v2
* elerch/S3 Filesystem

The Java SDK version two carries with it a significant amount of additional
dependencies. These had a significant effect on the ability of Cromwell to
compile, but ultimately the overlapping dependencies came down to:

* com.fasterxml.jackson.core (jackson-annotations):
    Jackson version was bumped to accomodate throughput
* org.slf4j (jcl-over-slf4j):
    Ignored by aws sdk in favor of the version bundled in Cromwell
* nettyHandler:
    Version bumped to accomodate AWS SDK
* sttpV:
    Version bumped to accomodate AWS SDK

While the AWS SDK v2 has nio capabilities, it does not include
a FileSystemProvider to allow Cromwell to operate seamlessly with AWS S3.
This was found later in the project, and as such, the filesystem provider
is not integrated with the Cromwell configuration system. The S3 FS was
forked from Udacity's provider, which was based on AWS SDK for Java version 1.
Of particular note is the fact that all API configuration for the provider
is performed through environment variables. As such, configuration is possible,
but currently disjoint from Cromwell proper.

Authentication
--------------

Authentication is required for both the Cromwell AWS Backend and the S3
Filesystem provider. By default, in both cases, the default credential provider
chain is followed. AWS tools all follow the same prioritized set of
checks for access key/secret key (and optionally token), with the exception
of #2 below.

[Default credential provider chain](https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html)
1. Configured permissions (environment)
2. (NOTE BELOW) Java properties
3. Access key/secret key as defined in $HOME/.aws/credentials and config
4. Container role
5. EC2 Instance role

NOTE: The Java properties check is specific and unique to the Java SDK. This
      does not apply to other SDKs or tools (e.g. the CLI).

Normally customers will be using an EC2 Instance role (recommended) or file configuration
as described in #3 (not recommended in production).

Permissions
-----------

Within AWS, everything must be authorized. This is a consistent rule, and as
such, AWS Services themselves are not immune to the rule. Therefore, customers
of AWS are responsible for granting services access to APIs within their account.
The flow described below represents the permissions needed by each stage, from 
Cromwell server through the task running. This includes the permissions needed for
the AWS Services involved in the processing of the work.

```text
+----------------------------+
|                            |  s3:GetObject on bucket for workflow and script bucket
|                            |  s3:ListObjects on script bucket
|                            |  s3:PutObject on script bucket
|          Cromwell          |  batch:RegisterTaskDefinition
|                            |  batch:SubmitJob
|                            |  batch:DescribeJobs
|                            |  batch:DescribeJobDefinitions
+-------------+--------------+
              |
              |
              |
+-------------v--------------+
|                            |  AWSBatchServiceRole managed policy - described at:
|          AWS Batch         |
|                            |     https://docs.aws.amazon.com/batch/latest/userguide/service_IAM_role.html
+-------------+--------------+
              |
              |
              |
+-------------v--------------+
|                            |  AWSServiceRoleForECS Service-linked role, documented at:
|                            |
| Elastic Container Service  |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/using-service-linked-roles.html
|                            |
|  (See discussion #1 below) |  AmazonEC2ContainerServiceAutoscaleRole managed policy - described at:
|                            |
|                            |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/autoscale_IAM_role.html
+-------------+--------------+
              |
              |
              |
+-------------v--------------+
|                            |
|                            |  AmazonEC2ContainerServiceforEC2Role managed policy, described at:
| ECS Agent (running on EC2) |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/instance_IAM_role.html (EC2)
|                            |    OR
|                            |  AmazonECSTaskExecutionRolePolicy managed policy, described at:  
|                            |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task_execution_IAM_role.html (Fargate)
+-------------+--------------+ 
              |
              |
              |
+-------------v--------------+
|                            |  Task Role permissions. These are user defined, but ecs-tasks.amazon.com must have sts:AssumeRole trust relationship defined. Documentation:
|       Task Container       |
|                            |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task_IAM_role.html
|                            |  s3:GetObject, s3:PutObject, s3:ListObjects
+----------------------------+
```


1. ECS has several sets of permissions for various items. AWS Batch, however,
   does not take advantage of certain features of ECS, most importantly
   ECS Services are out of scope of AWS Batch. ECS services require things
   like load balancing registration and DNS updates. While there is
   documentation regarding roles related to ECS services, these are
   irrelevant to the Cromwell use case.
2. Other than access to the main Cromwell bucket, the task container itself 
   does not need additional permissions unless
   the task in the WDL has been defined with a command that interfaces
   with AWS directly. This may include access to additional s3 buckets containing
   things like reference genome files.

   Task container permissions are currently supported through ECS and AWS
   Batch, but there is no configuration currently wired for the Cromwell
   AWS Backend to pass these settings to AWS. As such, the task container
   permissions must be managed by attaching a role to the EC2 Instance
   with permissions necessary for both the ECS Agent and the task container.

NOTE: ECS Agent permissions currently must use the permissions as outlined
      in the AmazonEC2ContainerServiceForEC2Role managed policy.

Future considerations
---------------------

AWS Batch Backend

* Should the 'disks' configuration be renamed to mount points or maybe
  ignored? This might make the wdl configuration more portable between backends
* The intent is to be able to override the queueArn between the
  default-runtime-attributes and the runtime attributes in the WDL. This
  appears broken at the moment.
* Job caching appears to be broken at the moment. Identical tasks need not be repeated
  if the results of a previous run of the task are still available.
* Retrying failed jobs is not currently attempted. Adding this would be beneficial 
  especially in conjunction with result caching
* Some S3 FS stuff can be removed: It looks like there’s a bunch of leftover
  unused code here from before having the nio implementation
  [here](https://github.com/broadinstitute/cromwell/tree/develop/filesystems/s3/src/main/scala/cromwell/filesystems/s3/batch)
  Also, I recommend adding another case for S3
  [here](https://github.com/broadinstitute/cromwell/blob/7a830a7fceaab9e7eaaa6802e58fe3cfd5d411a8/engine/src/main/scala/cromwell/engine/io/nio/NioFlow.scala#L95)
  otherwise files will be pulled down to be hashed instead of looking up the
  hash from S3 metadata. You might need to add the s3 Cromwell fs as a dependency
  to the engine to be able to do that.
* S3 Filesystem should be an official AWS thing (e.g. upplication or awslabs account)
* 8k container overrides limit should be bigger (ECS)
* We should understand why AWS Batch is putting everything in container overrides
* Authentication configuration consistency between S3FS and AWS Batch backend
* Full configuration of jobs in AWS Batch

Cromwell

* The style of integration with backends requires a lot of boilerplate code
  that I believe can be significantly reduced via heavier use of traits and
  supporting libraries
* There is a lot of data available within the backend if you know where
  to look, but if you don't, there is a lot of looking at inherited classes
  etc. to try to find them. Often, lack of code comments, etc and very similarly
  named variables can be confusing for folks confronting the code base for the
  first time. Workflow root vs call root is a great example. Adding additional
  comments and potentially context objects to encapsulate state would make
  the codebase more approachable.
* There is a significant amount of dependent libraries (with more added by
  the introduction of the S3 filesystem and the AWS SDK). Dependency management
  is challenging as a result. Adding significant new functionality is relatively
  painful when new dependencies are needed.
