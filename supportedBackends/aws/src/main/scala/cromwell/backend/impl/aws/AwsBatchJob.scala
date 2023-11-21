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

import cats.data.ReaderT._
import cats.data.{Kleisli, ReaderT}
import cats.effect.{Async, Timer}
import cats.syntax.all._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.impl.aws.io.AwsBatchWorkingDisk
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
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.{GetObjectRequest, HeadObjectRequest, NoSuchKeyException, PutObjectRequest}
import wdl4s.parser.MemoryUnit

import java.security.MessageDigest
import scala.concurrent.duration._
import scala.jdk.CollectionConverters._
import scala.util.{Random, Try}

/**
  *  The actual job for submission in AWS batch. `AwsBatchJob` is the primary interface to AWS Batch. It creates the
  *  necessary `AwsBatchJobDefinition`, then submits the job using the SubmitJob API.
  *
  *  @constructor Create an `AwsBatchJob` object capable of submitting, aborting and monitoring itself
  *  @param jobDescriptor the `BackendJobDescriptor` that is passed to the BackendWorkflowActor
  *  @param runtimeAttributes runtime attributes class (which subsequently pulls from config)
  *  @param commandLine command line to be passed to the job
  *  @param commandScript the `commandLine` with additional commands to setup the context and localize/ de-localize files
  */
final case class AwsBatchJob(jobDescriptor: BackendJobDescriptor, // WDL/CWL
                             runtimeAttributes: AwsBatchRuntimeAttributes, // config or WDL/CWL
                             commandLine: String, // WDL/CWL
                             commandScript: String, // WDL/CWL
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

  // values for container environment
  val AWS_MAX_ATTEMPTS: String = "AWS_MAX_ATTEMPTS"
  val AWS_MAX_ATTEMPTS_DEFAULT_VALUE: String = "14"
  val AWS_RETRY_MODE: String = "AWS_RETRY_MODE"
  val AWS_RETRY_MODE_DEFAULT_VALUE: String = "adaptive"

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

  /**
    * The goal of the reconfigured script is to do 3 things:
    * 1. make a directory "working_dir" and replace all "/cromwell_root" with that
    * 2. before the command line, insert the input copy command to stage input files
    * 3. at the end of the script sync all content with the s3 bucket
    */
  lazy val reconfiguredScript: String = {
    //this is the location of the aws cli mounted into the container by the ec2 launch template
    val s3Cmd = "/usr/local/aws-cli/v2/current/bin/aws s3"
    //internal to the container, therefore not mounted
    val workDir = "/tmp/scratch"
    //working in a mount will cause collisions in long running workers
    val replaced = commandScript.replace(AwsBatchWorkingDisk.MountPoint.pathAsString, workDir)
    val insertionPoint = replaced.indexOf("\n", replaced.indexOf("#!")) +1 //just after the new line after the shebang!

    /* generate a series of s3 copy statements to copy any s3 files into the container. We randomize the order
       so that large scatters don't all attempt to copy the same thing at the same time. */
    val inputCopyCommand = Random.shuffle(inputs.map {
      case input: AwsBatchFileInput if input.s3key.startsWith("s3://") && input.s3key.endsWith(".tmp") =>
        //we are localizing a tmp file which may contain workdirectory paths that need to be reconfigured
        s"""
           |$s3Cmd cp --no-progress ${input.s3key} $workDir/${input.local}
           |sed -i 's#${AwsBatchWorkingDisk.MountPoint.pathAsString}#$workDir#g' $workDir/${input.local}
           |""".stripMargin


      case input: AwsBatchFileInput if input.s3key.startsWith("s3://") =>
        s"$s3Cmd cp --no-progress ${input.s3key} ${input.mount.mountPoint.pathAsString}/${input.local}"
          .replace(AwsBatchWorkingDisk.MountPoint.pathAsString, workDir)

      case input: AwsBatchFileInput =>
        //here we don't need a copy command but the centaurTests expect us to verify the existence of the file
        val filePath = s"${input.mount.mountPoint.pathAsString}/${input.local.pathAsString}"
          .replace(AwsBatchWorkingDisk.MountPoint.pathAsString, workDir)

        s"test -e $filePath || echo 'input file: $filePath does not exist' && exit 1"

      case _ => ""
    }.toList).mkString("\n")

    // this goes at the start of the script after the #!
    val preamble =
      s"""
         |{
         |set -e
         |echo '*** LOCALIZING INPUTS ***'
         |if [ ! -d $workDir ]; then mkdir $workDir && chmod 777 $workDir; fi
         |cd $workDir
         |$inputCopyCommand
         |echo '*** COMPLETED LOCALIZATION ***'
         |set +e
         |}
         |""".stripMargin

    // the paths of the stdOut and stdErr
    val stdOut = dockerStdout.replace("/cromwell_root", workDir)
    val stdErr = dockerStderr.replace("/cromwell_root", workDir)

    //generate a series of s3 commands to delocalize artifacts from the container to storage at the end of the task
    val outputCopyCommand = outputs.map {
      case output: AwsBatchFileOutput if output.local.pathAsString.contains("*") => "" //filter out globs
      case output: AwsBatchFileOutput if output.name.endsWith(".list") && output.name.contains("glob-") =>

        val s3GlobOutDirectory = output.s3key.replace(".list", "")
        val globDirectory = output.name.replace(".list", "")
        /*
         * Need to process this list and de-localize each file if the list file actually exists
         * if it doesn't exist then 'touch' it so that it can be copied otherwise later steps will get upset
         * about the missing file
         */
        s"""
           |touch ${output.name}
           |$s3Cmd cp --no-progress ${output.name} ${output.s3key}
           |if [ -e $globDirectory ]; then $s3Cmd cp --no-progress $globDirectory $s3GlobOutDirectory --recursive --exclude "cromwell_glob_control_file"; fi
           |""".stripMargin

      case output: AwsBatchFileOutput if output.s3key.startsWith("s3://") && output.mount.mountPoint.pathAsString == AwsBatchWorkingDisk.MountPoint.pathAsString =>
        //output is on working disk mount
        s"""
           |$s3Cmd cp --no-progress $workDir/${output.local.pathAsString} ${output.s3key}
           |""".stripMargin
      case output: AwsBatchFileOutput =>
        //output on a different mount
        s"$s3Cmd cp --no-progress ${output.mount.mountPoint.pathAsString}/${output.local.pathAsString} ${output.s3key}"
      case _ => ""
    }.mkString("\n") + "\n" +
      s"""
         |if [ -f $workDir/${jobPaths.returnCodeFilename} ]; then $s3Cmd cp --no-progress $workDir/${jobPaths.returnCodeFilename} ${jobPaths.callRoot.pathAsString}/${jobPaths.returnCodeFilename} ; fi\n
         |if [ -f $stdErr ]; then $s3Cmd cp --no-progress $stdErr ${jobPaths.standardPaths.error.pathAsString}; fi
         |if [ -f $stdOut ]; then $s3Cmd cp --no-progress $stdOut ${jobPaths.standardPaths.output.pathAsString}; fi
         |""".stripMargin


    //insert the preamble at the insertion point and the postscript copy command at the end
    replaced.patch(insertionPoint, preamble, 0) +
      s"""
         |{
         |set -e
         |echo '*** DELOCALIZING OUTPUTS ***'
         |$outputCopyCommand
         |echo '*** COMPLETED DELOCALIZATION ***'
         |}
         |""".stripMargin
  }
  private def batch_file_s3_url(scriptBucketName: String, scriptKeyPrefix: String, scriptKey: String): String  = runtimeAttributes.fileSystem match {
       case AWSBatchStorageSystems.s3  => s"s3://${runtimeAttributes.scriptS3BucketName}/$scriptKeyPrefix$scriptKey"
       case _ => ""
  }

  private def generateEnvironmentKVPairs(scriptBucketName: String, scriptKeyPrefix: String, scriptKey: String): List[KeyValuePair] = {
    List(buildKVPair(AWS_MAX_ATTEMPTS, AWS_MAX_ATTEMPTS_DEFAULT_VALUE),
      buildKVPair(AWS_RETRY_MODE, AWS_RETRY_MODE_DEFAULT_VALUE),
      buildKVPair("BATCH_FILE_TYPE", "script"),
      buildKVPair("BATCH_FILE_S3_URL",batch_file_s3_url(scriptBucketName,scriptKeyPrefix,scriptKey)))
  }

  def submitJob[F[_]]()( implicit timer: Timer[F], async: Async[F]): Aws[F, SubmitJobResponse] = {

    val taskId = jobDescriptor.key.call.fullyQualifiedName + "-" + jobDescriptor.key.index + "-" + jobDescriptor.key.attempt

    //find or create the script in s3 to execute for s3 fileSystem
    val scriptKey =  runtimeAttributes.fileSystem match {
       case AWSBatchStorageSystems.s3  =>  findOrCreateS3Script(reconfiguredScript, runtimeAttributes.scriptS3BucketName)
       case _ => ""
    }

    if(runtimeAttributes.fileSystem == AWSBatchStorageSystems.s3) {
       val regex = "s3://([^/]*)/(.*)".r
       val regex(bucketName, key) = jobPaths.callExecutionRoot.toString
       writeReconfiguredScriptForAudit(reconfiguredScript, bucketName, key+"/reconfigured-script.sh")
    }


    val batch_script = runtimeAttributes.fileSystem match {
       case AWSBatchStorageSystems.s3  => s"s3://${runtimeAttributes.scriptS3BucketName}/$scriptKeyPrefix$scriptKey"
       case _  => commandScript
    }

    //calls the client to submit the job
    def callClient(definitionArn: String, awsBatchAttributes: AwsBatchAttributes): Aws[F, SubmitJobResponse] = {

      Log.info(s"Submitting taskId: $taskId, job definition : $definitionArn, script: $batch_script")

      val submit: F[SubmitJobResponse] =
        async.delay(batchClient.submitJob(
          SubmitJobRequest.builder()
            .jobName(sanitize(jobDescriptor.taskCall.fullyQualifiedName))
            .parameters(parameters.collect({ case i: AwsBatchInput => i.toStringString }).toMap.asJava)
            //provide job environment variables, vcpu and memory
            .containerOverrides(
              ContainerOverrides.builder
                .environment(

                  generateEnvironmentKVPairs(runtimeAttributes.scriptS3BucketName, scriptKeyPrefix, scriptKey): _*
                )
                .resourceRequirements(
                  ResourceRequirement.builder()
                    .`type`(ResourceType.VCPU)
                    .value(runtimeAttributes.cpu.value.toString)
                    .build(),
                  ResourceRequirement.builder()
                    .`type`(ResourceType.MEMORY)
                    .value(runtimeAttributes.memory.to(MemoryUnit.MB).amount.toInt.toString)
                    .build(),
                )
                .build()
            )
            .tags(runtimeAttributes.resourceTags.asJava)
            .jobQueue(runtimeAttributes.queueArn)
            .jobDefinition(definitionArn)
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

    Log.debug(s"s3 object name for script is calculated to be s3://$bucketName/$scriptKeyPrefix$key")

    try { //try and get the object

      s3Client.getObject( GetObjectRequest.builder().bucket(bucketName).key( scriptKeyPrefix + key).build )
      s3Client.headObject(HeadObjectRequest.builder()
        .bucket(bucketName)
        .key( scriptKeyPrefix + key )
        .build
      ).eTag().equals(key)

      // if there's no exception then the script already exists
      Log.debug(s"""Found script $bucketName/$scriptKeyPrefix$key""")
    } catch {
      case _: NoSuchKeyException =>  //this happens if there is no object with that key in the bucket
        val putRequest = PutObjectRequest.builder()
          .bucket(bucketName) //remove the "s3://" prefix
          .key(scriptKeyPrefix + key)
          .build

        s3Client.putObject(putRequest, RequestBody.fromString(commandLine))

        Log.debug(s"Created script $key")
    }
    key
  }

  private def writeReconfiguredScriptForAudit( reconfiguredScript: String, bucketName: String, key: String) = {
    val putObjectRequest = PutObjectRequest.builder().bucket(bucketName).key(key).build()
    s3Client.putObject(putObjectRequest, RequestBody.fromString(reconfiguredScript))
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
        case _ => commandScript
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

      Log.debug(s"Checking for existence of job definition called: $jobDefinitionName")

      val describeJobDefinitionRequest = DescribeJobDefinitionsRequest.builder()
        .jobDefinitionName( jobDefinitionName )
        .status("ACTIVE")
        .build()

      val describeJobDefinitionResponse = batchClient.describeJobDefinitions(describeJobDefinitionRequest)

      if ( !describeJobDefinitionResponse.jobDefinitions.isEmpty ) {
        //sort the definitions so that the latest revision is at the head
        val definitions = describeJobDefinitionResponse.jobDefinitions().asScala.toList.sortWith(_.revision > _.revision)

        //return the arn of the job
        definitions.head.jobDefinitionArn()
      } else {
        Log.debug(s"No job definition found. Creating job definition: $jobDefinitionName")

        // See:
        //
        // http://aws-java-sdk-javadoc.s3-website-us-west-2.amazonaws.com/latest/software/amazon/awssdk/services/batch/model/RegisterJobDefinitionRequest.Builder.html
        val definitionRequest = RegisterJobDefinitionRequest.builder
          .containerProperties(jobDefinition.containerProperties)
          .jobDefinitionName(jobDefinitionName)
          // See https://stackoverflow.com/questions/24349517/scala-method-named-type
          .`type`(JobDefinitionType.CONTAINER)
          .build

        Log.debug(s"Submitting definition request: $definitionRequest")

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
    runStatus <- RunStatus.fromJobStatus(statusString, jobId)
  } yield runStatus

  def detail(jobId: String): JobDetail = {
    //TODO: This client call should be wrapped in a cats Effect
    val describeJobsResponse = batchClient.describeJobs(DescribeJobsRequest.builder.jobs(jobId).build)

    val jobDetail = describeJobsResponse.jobs.asScala.headOption.
      getOrElse(throw new RuntimeException(s"Expected a job Detail to be present from this request: $describeJobsResponse and this response: $describeJobsResponse "))

    jobDetail
  }

  def rc(detail: JobDetail): Integer = {
    detail.container.exitCode
  }

  def output(detail: JobDetail): String = {
    val events: Seq[OutputLogEvent] = cloudWatchLogsClient.getLogEvents(GetLogEventsRequest.builder
      // http://aws-java-sdk-javadoc.s3-website-us-west-2.amazonaws.com/latest/software/amazon/awssdk/services/batch/model/ContainerDetail.html#logStreamName--
      .logGroupName(runtimeAttributes.logsGroup)
      .logStreamName(detail.container.logStreamName)
      .startFromHead(true)
      .build).events.asScala.toList
    val eventMessages = for ( event <- events ) yield event.message
    eventMessages mkString "\n"
  }

  def abort(jobId: String): TerminateJobResponse = {
    /*
     * Using Terminate here because it will work on jobs at any stage of their lifecycle whereas cancel will only work
     * on jobs that are not yet at the STARTING or RUNNING phase
     */
    batchClient.terminateJob(TerminateJobRequest.builder.jobId(jobId).reason("cromwell abort called").build())
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
      .append("commandScript", commandScript)
      .append("dockerRc", dockerRc).append("dockerStderr", dockerStderr).append("dockerStdout", dockerStdout)
      .append("inputs", inputs)
      .append("outputs", outputs)
      .append("jobPaths", jobPaths)
      .append("configRegion", configRegion)
      .append("awsAuthMode", optAwsAuthMode)
      .build
  }
}
