AWS Batch Backend Architecture
==============================

Overview
--------

The architecture of the code base follows very closely to the Google version.
Probably a little too closely, and lots of code was lifted from the Google
backend originally, then modified to work with AWS.

Fundamentally, Google Pipelines (a.k.a. JES) works pretty differently from
AWS Batch. In Pipelines, all the infrastructure is completely managed by Google,
while AWS Batch exposes that infrastructure to a large degree so that customers
can fine tune it as necessary.

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

AWS Batch fundamentally is a service to allow batch jobs to run easily and
efficiently. To use it effectively, however, you need to understand its own
technical stack. To create a job, you need a "Job Queue". That job queue allows
jobs to be scheduled onto one of several "Compte Environments". This can
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

ECS will run a special Amazon Machine Image (AMI) with a modified ecs agent,
a process for automatically expanding the Elastic Block Store (EBS) volumes,
and a proxy container. This will be discussed in more detail later.

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

1. The AwsBatchBackendLifecycleActorFactory class is configured by the user
   as the Cromwell backend. This factory provides an object from the
   AwsBatchAsyncBackendJobExecutionActor class to create and manage the job.
2. The AwsBatchAsyncBackendJobExecutionActor creates and manages the job.
   The job itself is encapsulated by the functionality in AwsBatchJob.
3. AwsBatchJob is the primary interface to AWS Batch. It creates the
   necessary AwsBatchJobDefinition, then submits the job using the SubmitJob
   API.
4. AwsBatchJobDefinition is responsible for the creation of the job definition.
   In AWS Batch, every job must have a definition. In the current backend,
   the job definition is a 1:1 with a job. Note that the job definition
   can be overridden by the SubmitJob, but in testing, container overrides
   have become problematic, possibly due to long environment variable values.
   Partly for this reason, all configuration information is stored in the
   Job Definition.

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
                ++-----------++
                 |           |
              Creates     Launches/
                 |        Monitors
                 |           |
+----------------v-+     +---v----------------+
|                  |     |                    |
|  Task Container  |     |   Proxy Container  |
|                  |     |                    |
+--------^---------+     +---------+----------+
         |                         |
         |                         |
         |                         |
         +---------Manages---------+
```

When a Cromwell task begins, the Cromwell backend will call the SubmitJob
API of AWS Batch. From there, the backend will call the AWS Batch DescribeJobs
API to provide status the the Cromwell engine as requested.

Once the job is Submitted in AWS Batch, one of the EC2 instances assigned
to the compute environment (a.k.a. ECS Cluster) with a running agent will
pick up the Cromwell Job/AWS Batch Job/ECS Task and run it. Importantly,
AWS Batch calls ECS' RunTask API when submitting the job. It uses the
task definition, and overrides both the command text and the environment
variables.

The agent is specifically modified for use with Cromwell. The agent will first
create the task container (but not start it), then duplicate all parameters
with a few exceptions (e.g. image name and task container ID environment
variable) as a proxy container. The agent, from then on, will launch and track
the proxy container rather than the task container.

Since the entire lifecycle of the task container is now the responsibility
of the proxy, the proxy can localize and delocalize files as it sees fit.
It uses 3 environment variables from the Task Container to do so:

* AWS_CROMWELL_OUTPUTS
* AWS_CROMWELL_INPUTS
* AWS_CROMWELL_INPUTS_GZ

These files describe the s3 key for the file, the location within the container,
and the mount point. The _GZ variant was added to avoid exceeding the AWS Batch
Limits, and is currently only used for inputs as fan-in is much more common
then fan out. However, this is easy to extend for outputs as well.

A few other variables are also important, and those are:

* AWS_CROMWELL_RC_FILE
* AWS_CROMWELL_STDOUT_FILE
* AWS_CROMWELL_STDERR_FILE

These files are placed in the correct location in S3 after task execution.

Logs for the task are added in the appropriate place so the normal assumptions
and console links operate as expected. In addition, the proxy is configured
by the agent to log side by side with the task container, such that if a
task is logging to the cloudwatch log stream abcdef-ghi, the proxy will
log to the log stream abcdef-ghi-proxy.

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

To address this issue and maximize the job size available, two changes were
made to the Cromwell backend. First, the AWS_CROMWELL_INPUTS environment
variable was gzip compressed and base64 encoded. Secondly, the command
itself was also compressed and encoded.

The AWS_CROMWELL_INPUTS environment variable describes all the inputs and
allows the proxy container to localize these files, copying them down from
S3. To compress this variable, code was added to the proxy to also look
for an AWS_CROMWELL_INPUTS_GZ variable. If found, the variable will be
base64 decoded and uncompressed, then processed as normal. The proxy is simply
a shell script, so this becomes a couple pipe commands input | base64 -d | zcat.

The command is much, much more difficult, but necessary as some instances of
MergeVCFs involve the command itself being upwards of 17KiB. As the container
command is processed entirely by the ecs-agent docker container, a second
change was made to the ecs-agent code base to inspect the command array. If
the magic string "gzipdata" was at element 0, the remaining array would be
copied to a new command array to be passed to docker. The last element of the
new array, however, would be base64 decoded, then uncompressed to the original.

At least in the Haplotype caller test, both the inputs and the command
text are highly compressible with lots of repeating text, so this scheme
should give a significant amount of headroom for other operations. The inputs
example can be easily extended by a small amount of change to the proxy and the
AwsBatchJobDefinition.scala code base. The downside is that from a Cromwell
user perspective, there is no specific limit on the size of commands as the
gzip operation will depend on the nuances of the underlying text.

In both cases (input and command text), this compression is conditional. The
idea here being that in the normal case, we want things in the console to
"look normal" to the customer. This also eases debugging. The downside, of
course, is that we add cyclomatic complexity and introduce the possibility
of things only failing when really big tasks are run.

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
    Jackson version was bumped to accomodate
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

Normally customers will be using an EC2 Instance role or file configuration
as described in #3.

Permissions
-----------

Within AWS, everything must be authorized. This is a consistent rule, and as
such, AWS Services themselves are not immune to the rule. Therefore, customers
of AWS are responsible for granting services access to APIs within their account.
The flow described below is represents from Cromwell through the task running,
and the permissions needed to run. This includes the permissions needed for
the AWS Services involved in the processing of the work.

```text
+----------------------------+
|                            |  s3:GetObject on bucket for workflow
|                            |  batch:RegisterTaskDefinition
|          Cromwell          |  batch:SubmitJob
|                            |  batch:DescribeJobs
|                            |  batch:DescribeJobDefinitions (retry only)
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
|                            |
+----------------------------+
```

Note that while the ECS Agent 
Additional permissions are needed for the proxy container, and that will be
discussed in #2 below.

1. ECS has several sets of permissions for various items. AWS Batch, however,
   does not take advantage of certain features of ECS, most importantly
   ECS Services are out of scope of AWS Batch. ECS services require things
   like load balancing registration and DNS updates. While there is
   documentation regarding roles related to ECS services, these are
   irrelevant to the Cromwell use case.
2. The task container itself does not need additional permissions unless
   the task in the WDL has been defined with a command that interfaces
   with AWS directly. However, the proxy container needs to send stdout,
   stderr, and the return code information to S3 for Cromwell to read after
   processing. Given that the task container and the proxy container currently
   share permissions, the task container must have s3:PutObject permissions
   for stdout/stderr/rc files for the workflow bucket.

   Task container permissions are currently supported through ECS and AWS
   Batch, but there is no configuration currently wired for the Cromwell
   AWS Backend to pass these settings to AWS. As such, the task container
   permissions must be managed by attaching a role to the EC2 Intance
   with permissions necessary for both the ECS Agent and the task container.

NOTE: ECS Agent permissions currently must use the permissions as outlined
      in the AmazonEC2ContainerServiceForEC2Role managed policy. Since
      Cromwell on AWS currently requires a custom AMI and modified ECS agent,
      Fargate is not an option and therefore EC2 permissions apply.

Future considerations
---------------------

AWS Batch Backend

* Should the 'disks' configuration be renamed to mount points or maybe
  ignored? This might make the wdl configuration more portable between backends
* The intent is to be able to override the queueArn between the
  default-runtime-attributes and the runtime attributes in the WDL. This
  appears broken at the moment.
* In retrospect it's unclear that the AWS Java SDK v2 was strictly needed,
  and because no S3 Filesystem provider existed, there was an impact to the
  speed of delivery. It did, however, position the backend better for the
  future.
* Job Definitions can be shared, but no attempt has been made to do so.
  Some hashing technique could be leveraged to re-use job definitions. However,
  in doing so there may be loss of information, as the only deterministic way
  to do this would be to make the job definition name a hash. This name
  is used in the Amazon Resource Name (ARN), and can thus be easily verified.
  Unfortunately, from a user perspective, the console would simply show a bunch
  of hashes rather than anything super useful.
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
