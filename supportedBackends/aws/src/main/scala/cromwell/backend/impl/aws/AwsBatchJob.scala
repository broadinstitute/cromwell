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

import cats.data.{Kleisli, ReaderT}
import cats.data.ReaderT._
import cats.implicits._
import cromwell.core.ExecutionIndex._

import scala.language.higherKinds
import cats.effect.{Async, Timer}
import com.amazonaws.services.s3.AmazonS3Client
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model._
import software.amazon.awssdk.services.cloudwatchlogs.CloudWatchLogsClient
import software.amazon.awssdk.services.cloudwatchlogs.model.{GetLogEventsRequest, OutputLogEvent}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.JobPaths
import cromwell.cloudsupport.aws.AwsConfiguration
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import org.slf4j.{Logger, LoggerFactory}
import fs2.Stream
import software.amazon.awssdk.core.sync.RequestBody
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.{HeadObjectRequest, ListObjectsV2Request, PutObjectRequest}

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.util.Try

/**
  *  The actual job for submission in AWS batch. `AwsBatchJob` is the primary interface to AWS Batch. It creates the
  *  necessary `AwsBatchJobDefinition`, then submits the job using the SubmitJob API.
  *
  *  Currently, each job will have its own job definition and queue. Support for separation and reuse of job
  *  definitions will be provided later.
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
  // TODO: Auth, endpoint
  lazy val client: BatchClient = {
    val builder = BatchClient.builder()
    optAwsAuthMode.foreach { awsAuthMode =>
      builder.credentialsProvider(awsAuthMode.provider())
    }
    configRegion.foreach(builder.region)
    builder.build
  }
  lazy val logsclient: CloudWatchLogsClient = {
    val builder = CloudWatchLogsClient.builder()
    configRegion.foreach(builder.region)
    builder.build
  }

  lazy val s3Client: S3Client = {
    val s3Client = AmazonS3Client.builder.build
  }

  lazy val reconfiguredScript: String = {
    // We'll use the MD5 of the dockerRc so the boundary is "random" but consistent
    val boundary = MessageDigest.getInstance("MD5").digest(dockerRc.getBytes).map("%02x".format(_)).mkString

    // NOTE: We are assuming the presence of a volume named "local-disk".
    //       This requires a custom AMI with the volume defined. But, since
    //       we need a custom AMI anyway for any real workflow, we just need
    //       to make sure this requirement is documented.
    //
    //       This MIME format is technically no longer necessary as the
    //       proxy docker container will manage the stdout/stderr/rc
    //       stuff being copied correctly, but this can still be useful
    //       later, so I'm leaving it here for potential future use.
    script.concat(s"""
    |echo "MIME-Version: 1.0
    |Content-Type: multipart/alternative; boundary="$boundary"
    |
    |--$boundary
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="rc.txt"
    |"
    |cat $dockerRc
    |echo "--$boundary
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stdout.txt"
    |"
    |cat $dockerStdout
    |echo "--$boundary
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stderr.txt"
    |"
    |cat $dockerStderr
    |echo "--$boundary--"
    |exit 0
    """).stripMargin
  }

  def submitJob[F[_]]()( implicit timer: Timer[F], async: Async[F]): Aws[F, SubmitJobResponse] = {

    val taskId = jobDescriptor.key.call.fullyQualifiedName + "-" + jobDescriptor.key.index + "-" + jobDescriptor.key.attempt
    val workflow = jobDescriptor.workflowDescriptor
    val uniquePath = workflow.callable.name + "/" +
                     jobDescriptor.taskCall.localName + "/" +
                     workflow.id + "/" +
                     jobDescriptor.key.index.fromIndex + "/" +
                     jobDescriptor.key.attempt
    Log.info(s"""Submitting job to AWS Batch""")
    Log.info(s"""dockerImage: ${runtimeAttributes.dockerImage}""")
    Log.info(s"""jobQueueArn: ${runtimeAttributes.queueArn}""")
    Log.info(s"""taskId: $taskId""")
    Log.info(s"""hostpath root: $uniquePath""")

    //find or create the script in s3 to execute
    val scriptKey = findOrCreateS3Script(commandLine, runtimeAttributes.scriptS3BucketName)

    def callClient(definitionArn: String, awsBatchAttributes: AwsBatchAttributes): Aws[F, SubmitJobResponse] = {
      val submit: F[SubmitJobResponse] =
        async.delay(client.submitJob(
          SubmitJobRequest.builder()
            .jobName(sanitize(jobDescriptor.taskCall.fullyQualifiedName))
            .parameters(parameters.collect({ case i: AwsBatchInput => i.toStringString }).toMap.asJava)
            .containerOverrides( ContainerOverrides.builder.environment(
                KeyValuePair.builder.name("BATCH_FILE_TYPE").value("script").build,
                KeyValuePair.builder.name("BATCH_FILE_S3_URL").value(s"""s3://${runtimeAttributes.scriptS3BucketName}/$scriptKey""").build
              ).build
            )
            .jobQueue(runtimeAttributes.queueArn)
            .jobDefinition(definitionArn).build
        ))

      ReaderT.liftF(
        Stream.retry(submit, 0.millis, duration => duration.plus(duration), awsBatchAttributes.submitAttempts.value, {
          // RegisterJobDefinition is eventually consistent, so it may not be there
          case e: ClientException => e.statusCode() == 404
          case _ => false
        }).compile.last.map(_.get)) //if successful there is guaranteed to be a value emitted, hence we can .get this option
    }

    (createDefinition[F](s"""${workflow.callable.name}-${jobDescriptor.taskCall.localName}""", uniquePath) product Kleisli.ask[F, AwsBatchAttributes]).
      flatMap((callClient _).tupled)
  }

  /**
    * is there an existing appropriate definition that I can use, if not create one
    */
  private def findOrRegisterJobDefinition(): Unit ={

    //is there an active definition of the same name?
    val jobDefinitionName = s"""cromwell_def_${runtimeAttributes.dockerImage}"""

    val jobDefinitionRequest = DescribeJobDefinitionsRequest.builder()
      .jobDefinitionName( sanitize(jobDefinitionName) )
      .status( "ACTIVE" )
      .build()

    val existingJobDefinitions = client.describeJobDefinitions(jobDefinitionRequest).jobDefinitions()

    if( existingJobDefinitions.isEmpty ) {
      //can't find one of the same name so make it
      Log.info(s"Creating $jobDefinitionName job definition")

      val jobRegistrationResult: RegisterJobDefinitionResponse = registerJobDefinition(jobDefinitionName)

      Log.info(s"Created job definition named ${jobRegistrationResult.jobDefinitionName()} " +
        s"with arn ${jobRegistrationResult.jobDefinitionArn}")

    } else {
      //get the latest revision
      val jobDef = existingJobDefinitions.asScala.toList.sortWith( _.revision > _.revision ).head

      println(s"Found job definition ${ jobDef.jobDefinitionName } revision ${ jobDef.revision}")

      //todo check that the definition that was found by name matches what would have been created in code??

      // the definition should have the same name and the latest version should have the expected:
      // 1. image
      // 2. command (fetch and run)
      // 3. volumes:
      //   "volumes": [
      //            {
      //                "host": {
      //                    "sourcePath": "/usr/local/bin/fetch_and_run.sh"
      //                },
      //                "name": "fetchAndRunScript"
      //            },
      //            {
      //                "host": {
      //                    "sourcePath": "/usr/local/aws-cli"
      //                },
      //                "name": "awsCliHome"
      //            }
      //   ]
      // 4. mount points:
      //   "mountPoints": [
      //            {
      //                "containerPath": "/var/scratch/fetch_and_run.sh",
      //                "readOnly": true,
      //                "sourceVolume": "fetchAndRunScript"
      //            },
      //            {
      //                "containerPath": "/usr/local/aws-cli",
      //                "readOnly": true,
      //                "sourceVolume": "awsCliHome"
      //            }
      //   ]

    }
  }

//  private def registerJobDefinition( jobDefinitionName: String ): RegisterJobDefinitionResponse = {
//    val registerJobDefinitionRequest = RegisterJobDefinitionRequest.builder()
//      .jobDefinitionName(jobDefinitionName)
//      // type is reserved in scala so java methods called type() need to be escaped
//      .`type`(JobDefinitionType.CONTAINER)
//      .containerProperties(ContainerProperties.builder()
//        //.jobRoleArn( jobRoleArn )
//        .image(runtimeAttributes.dockerImage)
//        // mount the fetch and run script. The script is made available by the launch template user data
//        .volumes(Volume.builder()
//          .name("fetchAndRunScript")
//          .host(Host.builder().sourcePath("/usr/local/bin/fetch_and_run.sh").build())
//          .build())
//        .mountPoints(MountPoint.builder()
//          .readOnly(true)
//          .sourceVolume("fetchAndRunScript")
//          .containerPath("/var/scratch/fetch_and_run.sh")
//          .build())
//        //mount the aws cli v2 distribution so the container can access it
//        .volumes(Volume.builder()
//          .name("awsCliHome")
//          .host(Host.builder().sourcePath("/usr/local/aws-cli").build())
//          .build())
//        .mountPoints(MountPoint.builder()
//          .readOnly(true)
//          .sourceVolume("awsCliHome")
//          .containerPath("/usr/local/aws-cli")
//          .build())
//        //the command the container should run
//        .command("/var/scratch/fetch_and_run.sh")
//        .build()
//      )
//      .build()
//
//    client.registerJobDefinition(registerJobDefinitionRequest)
//  }


  /**
    * Performs an md5 digest the command line, checks in s3 bucket for that script, if it's not already there then persists it.
    *
    * @param commandLine the command line script to be executed
    * @param scriptS3BucketName the bucket that stores the scripts
    */
  private def findOrCreateS3Script(commandLine :String, scriptS3BucketName: String) :String = {

    val s3Client = S3Client.create()

    val key = MessageDigest.getInstance("MD5")
      .digest(commandLine.getBytes())
      .foldLeft("")(_ + "%02x".format(_))

    val bucketName = scriptS3BucketName;

    if (s3Client.headObject(HeadObjectRequest.builder()
      .bucket(bucketName)
      .ifMatch(key)
      .build()
    ).eTag().equals(key)) {
      //the script already exists
      Log.info(s"""Found script s3://$bucketName/$key""")

    } else {

      Log.info(s"Script $key not found in bucket $bucketName. Creating script with content:\n$commandLine")
      s3Client.putObject(PutObjectRequest.builder()
        .bucket(bucketName)
        .key(key)
        .contentMD5(key)
        .build(),
        RequestBody.fromString(commandLine))

      Log.info(s"Created script $key")

    }

    key
  }

  /** Creates a job definition in AWS Batch
   *
   *  @param name Job definition name
   *  @return Arn for newly created job definition
   *
   */
  private def createDefinition[F[_]](name: String,
                                     taskId: String)(
                                     implicit async: Async[F],
                                     timer: Timer[F]): Aws[F, String] = ReaderT { awsBatchAttributes =>
    //todo check if there is already a suitable definition

    val jobDefinitionBuilder = StandardAwsBatchJobDefinitionBuilder
    val commandStr = awsBatchAttributes.fileSystem match {
      case AWSBatchStorageSystems.s3  => reconfiguredScript
      case _ => script
    }
    val jobDefinitionContext = AwsBatchJobDefinitionContext(runtimeAttributes,
                                                            taskId,
                                                            commandStr,
                                                            dockerRc,
                                                            dockerStdout,
                                                            dockerStderr,
                                                            jobDescriptor,
                                                            jobPaths,
                                                            inputs,
                                                            outputs)

    Log.info(s"Creating job definition for task: ${taskId} with command: ${commandStr}")

    val jobDefinition = jobDefinitionBuilder.build(jobDefinitionContext)

    // See:
    //
    // http://aws-java-sdk-javadoc.s3-website-us-west-2.amazonaws.com/latest/software/amazon/awssdk/services/batch/model/RegisterJobDefinitionRequest.Builder.html
    val definitionRequest = RegisterJobDefinitionRequest.builder
                              .containerProperties(jobDefinition.containerProperties)
                              .jobDefinitionName(sanitize(name))
                              // See https://stackoverflow.com/questions/24349517/scala-method-named-type
                              .`type`(JobDefinitionType.CONTAINER)
                              .build

    Log.info(s"Submitting definition request: ${definitionRequest}")
    val submit = async.delay(client.registerJobDefinition(definitionRequest).jobDefinitionArn)

    //Log.info(s"Submitted: ${submit.show}")

    val retry: F[String] = Stream.retry(submit, 0.millis, _ * 2, awsBatchAttributes.createDefinitionAttempts.value, {
      // RegisterJobDefinition throws 404s every once in a while
      case e: ClientException => e.statusCode() == 404
      case _ => false
    }).compile.last.map(_.get)

    async.recoverWith(retry) {
      case e: ClientException if e.statusCode() == 409 =>
        // This could be a problem here, as we might have registerJobDefinition
        // but not describeJobDefinitions permissions. We've changed this to a
        // warning as a result of this potential
        Log.warn("Job definition already exists. Performing describe and retrieving latest revision.")
        async.
          delay(client.describeJobDefinitions(DescribeJobDefinitionsRequest.builder().jobDefinitionName(sanitize(name)).build())).
          map(_.jobDefinitions().asScala).
          map(definitions => {
            Log.info(s"Latest job definition revision is: ${definitions.last.jobDefinitionName()} with arn: ${definitions.last.jobDefinitionArn()}")
            definitions.last.jobDefinitionArn()})
    }
  }

  /** Sanitizes a job and job definition name
    *
    *  @param name Job or Job definition name
    *  @return Sanitized name
    *
    */
  private def sanitize(name: String): String =
    // Up to 128 letters (uppercase and lowercase), numbers, hyphens, and underscores are allowed.
    // We'll replace all invalid characters with an underscore
    name.replaceAll("[^A-Za-z0-9_\\-]", "_").slice(0,128)


  /** Gets the status of a job by its Id, converted to a RunStatus
   *
   *  @param jobId Job ID as defined in AWS Batch
   *  @return Current RunStatus
   *
   */
  def status(jobId: String): Try[RunStatus] = for {
    statusString <- Try(detail(jobId).status)
    runStatus <- RunStatus.fromJobStatus(statusString, jobId)
  } yield runStatus

  def detail(jobId: String): JobDetail = {
    //TODO: This client call should be wrapped in a cats Effect
     val describeJobsResponse = client.describeJobs(DescribeJobsRequest.builder.jobs(jobId).build)

     describeJobsResponse.jobs.asScala.headOption.
       getOrElse(throw new RuntimeException(s"Expected a job Detail to be present from this request: $describeJobsResponse and this response: $describeJobsResponse "))
  }

  //TODO: unused at present
  def rc(detail: JobDetail): Integer = {
     detail.container.exitCode
  }

  def output(detail: JobDetail): String = {
     val events: Seq[OutputLogEvent] = logsclient.getLogEvents(GetLogEventsRequest.builder
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
    client.cancelJob(CancelJobRequest.builder.jobId(jobId).reason("cromwell abort called").build)
  }
}
