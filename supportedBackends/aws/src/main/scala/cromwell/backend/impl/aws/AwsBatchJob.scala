/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
package cromwell.backend.impl.aws

import java.security.MessageDigest

import cats.data.ReaderT._
import cats.data.{Kleisli, ReaderT}
import cats.effect.{Async, Timer}
import cats.implicits._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.JobPaths
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import fs2.Stream
import org.apache.commons.lang3.builder.{ToStringBuilder, ToStringStyle}
import org.slf4j.{Logger, LoggerFactory}
import software.amazon.awssdk.core.sync.RequestBody
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model._
import software.amazon.awssdk.services.cloudwatchlogs.CloudWatchLogsClient
import software.amazon.awssdk.services.cloudwatchlogs.model.{GetLogEventsRequest, OutputLogEvent}
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.ecs.model.DescribeContainerInstancesRequest
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.{GetObjectRequest, HeadObjectRequest, NoSuchKeyException, PutObjectRequest}
import wdl4s.parser.MemoryUnit

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.language.higherKinds
import scala.util.Try

/**
  *  The actual job for submission in AWS batch. `AwsBatchJob` is the primary interface to AWS Batch. It creates the
  *  necessary `AwsBatchJobDefinition`, then submits the job using the SubmitJob API.
  *
  *  @constructor Create an `AwsBatchJob` object capable of submitting, aborting and monitoring itself
  *  @param jobDescriptor the `BackendJobDescriptor` that is passed to the BackendWorkflowActor
  *  @param runtimeAttributes runtime attributes class (which subsequently pulls from config)
  *  @param commandLine command line to be passed to the job
  */
final case class AwsBatchJob(jobDescriptor: BackendJobDescriptor, // WDL/CWL
                             runtimeAttributes: AwsBatchRuntimeAttributes, // config or WDL/CWL
                             commandLine: String, // WDL/CWL
                             script: String, // WDL/CWL
                             dockerRc: String, // Calculated from StandardAsyncExecutionActor
                             dockerStdout: String, // Calculated from StandardAsyncExecutionActor
                             dockerStderr: String, // Calculated from StandardAsyncExecutionActor
                             inputs: Set[AwsBatchInput],
                             outputs: Set[AwsBatchFileOutput],
                             jobPaths: JobPaths, // Based on config, calculated in Job Paths, key to all things outside container
                             parameters: Seq[AwsBatchParameter],
                             configRegion: Option[Region],
                             optAwsAuthMode: Option[AwsAuthMode] = None
                             ) {

  val Log: Logger = LoggerFactory.getLogger(AwsBatchJob.getClass)

  //this will be the "folder" that scripts will live in (underneath the script bucket)
  val scriptKeyPrefix = "scripts/"

  // TODO: Auth, endpoint
  lazy val batchClient: BatchClient = {
    val builder = BatchClient.builder()
    configureClient(builder, optAwsAuthMode, configRegion)
  }
  lazy val cloudWatchLogsClient: CloudWatchLogsClient = {
    val builder = CloudWatchLogsClient.builder()
    configureClient(builder, optAwsAuthMode, configRegion)
  }

  lazy val s3Client: S3Client = {
    val builder = S3Client.builder()
    configureClient(builder, optAwsAuthMode, configRegion)
  }


  lazy val reconfiguredScript: String = {
    //this is the location of the aws cli mounted into the container by the ec2 launch template
    val s3Cmd = "/usr/local/aws-cli/v2/current/bin/aws s3"
    //val dockerRootDir = runtimeAttributes.disks.map(_.mountPoint.toString).head
    val dockerRootDir = "."

    //generate a series of s3 copy statements to copy any s3 files into the container
    val inputCopyCommand = inputs.map {
      case input: AwsBatchFileInput if input.s3key.startsWith("s3://") => s"$s3Cmd cp ${input.s3key} $dockerRootDir/${input.local}"
      case _ => ""
    }.mkString("\n|") //the pipe is needed for the margin (which gets stripped out)

    val outputCopyCommand = outputs.map {
      case output: AwsBatchFileOutput if output.s3key.startsWith("s3://") => s"$s3Cmd cp ${output.local} ${output.s3key}"
      case _ => ""
    }.mkString("\n|")

    //this is a temporary work around until we can surgically remove the "local-disk" volume
    val replacementCommandLine = commandLine.replaceAllLiterally("/cromwell_root", dockerRootDir)

    s"""
    |#! /bin/sh
    |
    |#exit on any line causing an error
    |set -e
    |
    |#trap any error and handle with the final() function
    |trap 'final $$? $$LINENO' ERR
    |
    |final() {
       #create or overwrite rc.txt, something failedf
    |  echo $$1 > rc.txt
    |  echo "Error $$1 occurred on line $$2"
    |}
    |
    |$inputCopyCommand
    |
    |touch stdout.log && touch stderr.log
    |
    |$replacementCommandLine | tee > stdout.log 2> stderr.log
    |return_code=$$?
    |echo $$return_code > rc.txt
    |
    |cat stdout.log && cat stderr.log >&2
    |
    |$outputCopyCommand
    |
    |$s3Cmd cp stdout.log $${AWS_CROMWELL_CALL_ROOT}/${jobPaths.defaultStdoutFilename}
    |$s3Cmd cp stderr.log $${AWS_CROMWELL_CALL_ROOT}/${jobPaths.defaultStderrFilename}
    |$s3Cmd cp rc.txt $${AWS_CROMWELL_CALL_ROOT}/${jobPaths.returnCodeFilename}
    |exit $$return_code
    |""".stripMargin
  }


  def submitJob[F[_]]()( implicit timer: Timer[F], async: Async[F]): Aws[F, SubmitJobResponse] = {

    val taskId = jobDescriptor.key.call.fullyQualifiedName + "-" + jobDescriptor.key.index + "-" + jobDescriptor.key.attempt

    //find or create the script in s3 to execute
    val scriptKey = findOrCreateS3Script(reconfiguredScript, runtimeAttributes.scriptS3BucketName)

    //calls the client to submit the job
    def callClient(definitionArn: String, awsBatchAttributes: AwsBatchAttributes): Aws[F, SubmitJobResponse] = {

      Log.info(s"""Submitting job to AWS Batch
                  |  dockerImage: ${runtimeAttributes.dockerImage}
                  |  jobQueueArn: ${runtimeAttributes.queueArn}
                  |  taskId: $taskId
                  |  job definition arn: $definitionArn
                  |  executionScript: $scriptKey
                  |  """.stripMargin)

      val outputinfo = outputs.map(o => "%s,%s,%s,%s".format(o.name, o.s3key, o.local, o.mount))
        .mkString(";")
      val inputinfo = inputs.collect{case i: AwsBatchFileInput => i}
        .map(i => "%s,%s,%s,%s".format(i.name, i.s3key, i.local, i.mount))
        .mkString(";")

      val submit: F[SubmitJobResponse] =
        async.delay(batchClient.submitJob(
          SubmitJobRequest.builder()
            .jobName(sanitize(jobDescriptor.taskCall.fullyQualifiedName))
            .parameters(parameters.collect({ case i: AwsBatchInput => i.toStringString }).toMap.asJava)

            //provide job environment variables, vcpu and memory
            .containerOverrides(
              ContainerOverrides.builder
                .environment(
                  buildKVPair("BATCH_FILE_TYPE", "script"),
                  buildKVPair("BATCH_FILE_S3_URL",
                    s"s3://${runtimeAttributes.scriptS3BucketName}/$scriptKeyPrefix$scriptKey"),
                  buildKVPair("AWS_CROMWELL_CALL_ROOT", jobPaths.callExecutionRoot.toString),
                  buildKVPair("AWS_CROMWELL_WORKFLOW_ROOT", jobPaths.workflowPaths.workflowRoot.toString),
                  gzipKeyValuePair("AWS_CROMWELL_INPUTS", inputinfo),
                  buildKVPair("AWS_CROMWELL_OUTPUTS", outputinfo),
                  buildKVPair("AWS_CROMWELL_STDOUT_FILE", dockerStdout),
                  buildKVPair("AWS_CROMWELL_STDERR_FILE", dockerStderr),
                  buildKVPair("AWS_CROMWELL_RC_FILE", dockerRc),
                )
                .memory(runtimeAttributes.memory.to(MemoryUnit.MB).amount.toInt)
                .vcpus(runtimeAttributes.cpu.##).build
            )
            .jobQueue(runtimeAttributes.queueArn)
            .jobDefinition(definitionArn)
            // todo add JobRetryStrategy with retry attempts
            .build
        ))

      ReaderT.liftF(
        Stream.retry(submit, 0.millis, duration => duration.plus(duration), awsBatchAttributes.submitAttempts.value, {
          // RegisterJobDefinition is eventually consistent, so it may not be there
          case e: ClientException => e.statusCode() == 404
          case _ => false
        }).compile.last.map(_.get)) //if successful there is guaranteed to be a value emitted, hence we can .get this option
    }

    (findOrCreateDefinition[F]() product Kleisli.ask[F, AwsBatchAttributes]).flatMap((callClient _).tupled)
  }


  /**
    * Performs an md5 digest the script, checks in s3 bucket for that script, if it's not already there then persists it.
    *
    * @param commandLine the command line script to be executed
    * @param scriptS3BucketName the bucket that stores the scripts
    * @return the name of the script that was found or created
    */
  private def findOrCreateS3Script(commandLine :String, scriptS3BucketName: String) :String = {

    val bucketName = scriptS3BucketName

    val key = MessageDigest.getInstance("MD5")
      .digest(commandLine.getBytes())
      .foldLeft("")(_ + "%02x".format(_))

    Log.info(s"s3 object name for script is calculated to be $bucketName/$scriptKeyPrefix$key")

    try { //try and get the object

      s3Client.getObject( GetObjectRequest.builder().bucket(bucketName).key( scriptKeyPrefix + key).build )
      s3Client.headObject(HeadObjectRequest.builder()
        .bucket(bucketName)
        .key( scriptKeyPrefix + key )
        .build
      ).eTag().equals(key)

      // if there's no exception then the script already exists
      Log.info(s"""Found script $bucketName/$scriptKeyPrefix$key""")
    } catch {
      case _: NoSuchKeyException =>  //this happens if there is no object with that key in the bucket
        Log.info(s"Script $key not found in bucket $bucketName. Creating script with content:\n$commandLine")

        val putRequest = PutObjectRequest.builder()
          .bucket(bucketName) //remove the "s3://" prefix
          .key(scriptKeyPrefix + key)
          .build

        s3Client.putObject(putRequest, RequestBody.fromString(commandLine))

        Log.info(s"Created script $key")
    }
    key
  }

  /** Creates a job definition in AWS Batch
   *
   * @return Arn for newly created job definition
   *
   */
  private def findOrCreateDefinition[F[_]]()
                                    (implicit async: Async[F], timer: Timer[F]): Aws[F, String] = ReaderT { awsBatchAttributes =>

    // this is a call back that is executed below by the async.recoverWithRetry(retry)
    val submit = async.delay({

      val commandStr = awsBatchAttributes.fileSystem match {
        case AWSBatchStorageSystems.s3 => reconfiguredScript
        case _ => script
      }
      val jobDefinitionContext = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = commandStr,
        dockerRcPath = dockerRc,
        dockerStdoutPath = dockerStdout,
        dockerStderrPath = dockerStderr,
        jobDescriptor = jobDescriptor,
        jobPaths = jobPaths,
        inputs = inputs,
        outputs = outputs)

      val jobDefinitionBuilder = StandardAwsBatchJobDefinitionBuilder
      val jobDefinition = jobDefinitionBuilder.build(jobDefinitionContext)


      //check if there is already a suitable definition based on the calculated job definition name
      val jobDefinitionName = jobDefinition.name

      Log.info(s"Checking for existence of job definition called: $jobDefinitionName")

      val describeJobDefinitionRequest = DescribeJobDefinitionsRequest.builder()
        .jobDefinitionName( jobDefinitionName )
        .status("ACTIVE")
        .build()

      val describeJobDefinitionResponse = batchClient.describeJobDefinitions(describeJobDefinitionRequest)

      if ( !describeJobDefinitionResponse.jobDefinitions.isEmpty ) {

        Log.info(s"Found job definition $jobDefinitionName, getting arn for latest version")

        //sort the definitions so that the latest revision is at the head
        val definitions = describeJobDefinitionResponse.jobDefinitions().asScala.toList.sortWith(_.revision > _.revision)

        Log.info(s"Latest job definition revision is: ${definitions.head.revision()} with arn: ${definitions.head.jobDefinitionArn()}")

        //return the arn of the job
        definitions.head.jobDefinitionArn()
      } else {

        //no definition found. create one
        Log.info(s"No job definition found")

        Log.info(s"Creating job definition: $jobDefinitionName")

        // See:
        //
        // http://aws-java-sdk-javadoc.s3-website-us-west-2.amazonaws.com/latest/software/amazon/awssdk/services/batch/model/RegisterJobDefinitionRequest.Builder.html
        val definitionRequest = RegisterJobDefinitionRequest.builder
          .containerProperties(jobDefinition.containerProperties)
          .jobDefinitionName(jobDefinitionName)
          // See https://stackoverflow.com/questions/24349517/scala-method-named-type
          .`type`(JobDefinitionType.CONTAINER)
          .build

        Log.info(s"Submitting definition request: $definitionRequest")

        val response: RegisterJobDefinitionResponse = batchClient.registerJobDefinition(definitionRequest)
        Log.info(s"Definition created: $response")
        response.jobDefinitionArn()
      }
    })


    // a function to retry submissions, returns a higher kind parameterized on a String (where the String is an arn)
    val retry: F[String] = Stream.retry(
      fo = submit, //the value to attempt to get
      delay = 0.millis, //how long to wait
      nextDelay = _ * 2, //how long to back off after a failure
      maxAttempts = awsBatchAttributes.createDefinitionAttempts.value, //how many times to try
      retriable = { //a function to say if we should retry or not
        // RegisterJobDefinition throws 404s every once in a while
        case e: ClientException => e.statusCode() == 404 || e.statusCode() == 409
        // a 409 means an eventual consistency collision has happened, most likely during a scatter.
        // Just wait and retry as job definition names are canonical and if another thread succeeds in making one then
        // that will be used and if there really isn't one, then the definition will be created.
        case _ => false  //don't retry other cases
      }
    ).compile.last.map(_.get)

    // attempt to register the job definition
    async.recoverWith(submit){
      case e: ClientException if e.statusCode == 404 ||
        e.statusCode == 409 || e.statusCode == 429  => retry  //probably worth trying again
    }
  }


  /** Gets the status of a job by its Id, converted to a RunStatus
   *
   *  @param jobId Job ID as defined in AWS Batch
   *  @return Current RunStatus
   *
   */
  def status(jobId: String): Try[RunStatus] = for {
    statusString <- Try(detail(jobId).status)
    batchJobContainerContext <- Try(batchJobContainerContext(jobId))
    _ <- Try(Log.info(s"Task ${jobDescriptor.key.call.fullyQualifiedName + "-" + jobDescriptor.key.index + "-" + jobDescriptor.key.attempt} in container context $batchJobContainerContext"))
    runStatus <- RunStatus.fromJobStatus(statusString, jobId)
  } yield runStatus

  def detail(jobId: String): JobDetail = {
    //TODO: This client call should be wrapped in a cats Effect
    val describeJobsResponse = batchClient.describeJobs(DescribeJobsRequest.builder.jobs(jobId).build)

    val jobDetail = describeJobsResponse.jobs.asScala.headOption.
      getOrElse(throw new RuntimeException(s"Expected a job Detail to be present from this request: $describeJobsResponse and this response: $describeJobsResponse "))

    jobDetail
  }

  /**
    * Return information about the container, ECS Cluster and EC2 instance that is (or was) hosting this job
    * @param jobId the id of the job for which you want the context
    * @return the context
    */
  def batchJobContainerContext(jobId: String): BatchJobContainerContext ={
    if (jobId == null) return BatchJobContainerContext("","",Seq.empty, Seq.empty)

    val containerInstanceArn = detail(jobId).container().containerInstanceArn()
    if(containerInstanceArn == null || containerInstanceArn.isEmpty) return BatchJobContainerContext(jobId,"",Seq.empty, Seq.empty)

    val describeJobQueuesResponse = batchClient.describeJobQueues( DescribeJobQueuesRequest.builder().jobQueues( runtimeAttributes.queueArn ).build())
    val computeEnvironments = describeJobQueuesResponse.jobQueues().asScala.head.computeEnvironmentOrder().asScala.map(_.computeEnvironment())
    val describeComputeEnvironmentsResponse = batchClient.describeComputeEnvironments( DescribeComputeEnvironmentsRequest.builder().computeEnvironments(computeEnvironments.asJava).build())
    val ecsClusterArns = describeComputeEnvironmentsResponse.computeEnvironments().asScala.map(_.ecsClusterArn())

    val ecsClient = configureClient(EcsClient.builder(), optAwsAuthMode, configRegion)

    val instanceIds: Seq[String] = ecsClusterArns.map(containerArn => ecsClient.describeContainerInstances(DescribeContainerInstancesRequest.builder().containerInstances(containerInstanceArn).cluster(containerArn).build()))
      .map(r => r.containerInstances().asScala).flatMap(_.map(_.ec2InstanceId()))

    BatchJobContainerContext(jobId, containerInstanceArn, ecsClusterArns, instanceIds)
  }

  case class BatchJobContainerContext(jobId: String, containerInstanceArn: String, ecsClusterArns: Seq[String], ec2InstanceIds: Seq[String]) {
    override def toString: String = {
      new ToStringBuilder(this, ToStringStyle.JSON_STYLE)
        .append("jobId", this.jobId)
        .append("containerInstanceArn", containerInstanceArn)
        .append("ecsClusterArns", ecsClusterArns)
        .append("ec2InstanceIds", ec2InstanceIds)
        .build()
    }
  }

  def rc(detail: JobDetail): Integer = {
     detail.container.exitCode
  }

  //todo: unused at present??
  def output(detail: JobDetail): String = {
     val events: Seq[OutputLogEvent] = cloudWatchLogsClient.getLogEvents(GetLogEventsRequest.builder
                                            // http://aws-java-sdk-javadoc.s3-website-us-west-2.amazonaws.com/latest/software/amazon/awssdk/services/batch/model/ContainerDetail.html#logStreamName--
                                            .logGroupName("/aws/batch/job")
                                            .logStreamName(detail.container.logStreamName)
                                            .startFromHead(true)
                                            .build).events.asScala
     val eventMessages = for ( event <- events ) yield event.message
     eventMessages mkString "\n"
  }

  //TODO: Wrap in cats Effect
  def abort(jobId: String): CancelJobResponse = {
    batchClient.cancelJob(CancelJobRequest.builder.jobId(jobId).reason("cromwell abort called").build)
  }

  /**
    * Generate a `String` describing the instance. Mainly for debugging
    * @return a description of the instance
    */
  override def toString: String = {
     new ToStringBuilder(this, ToStringStyle.JSON_STYLE)
      .append("jobDescriptor", jobDescriptor)
      .append("runtimeAttributes", runtimeAttributes)
      .append("commandLine", commandLine)
      .append("script", script)
      .append("dockerRc", dockerRc).append("dockerStderr", dockerStderr).append("dockerStdout", dockerStdout)
      .append("inputs", inputs)
      .append("outputs", outputs)
      .append("jobPaths", jobPaths)
      .append("configRegion", configRegion)
      .append("awsAuthMode", optAwsAuthMode)
      .build
  }
}
