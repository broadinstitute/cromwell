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
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model.
                                        {
                                          CancelJobRequest,
                                          CancelJobResponse,
                                          ClientException,
                                          DescribeJobDefinitionsRequest,
                                          DescribeJobsRequest,
                                          RegisterJobDefinitionRequest,
                                          SubmitJobRequest,
                                          SubmitJobResponse,
                                          JobDefinitionType,
                                          JobDetail
                                        }
import software.amazon.awssdk.services.cloudwatchlogs.CloudWatchLogsClient
import software.amazon.awssdk.services.cloudwatchlogs.model.GetLogEventsRequest
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.JobPaths
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.util.Try

object AwsBatchJob {
}

/** The actual job for submission in AWS batch. Currently, each job will
 *  have its own job definition and queue. Support for separation and reuse of job
 *  definitions will be provided later.
 *
 *  @constructor Create an AwsBatchJob object capable of submitting, aborting and monitoring itself
 *  @param jobDescriptor create a new person with a name and age.
 *  @param runtimeAttributes runtime attributes class (which subsequently pulls from config)
 *  @param commandLine command line to be passed to the job
 */
final case class AwsBatchJob(jobDescriptor: BackendJobDescriptor,           // WDL/CWL
                             runtimeAttributes: AwsBatchRuntimeAttributes,  // config or WDL/CWL
                             commandLine: String,                           // WDL/CWL
                             script: String,                                // WDL/CWL
                             dockerRc: String,                              // Calculated from StandardAsyncExecutionActor
                             dockerStdout: String,                          // Calculated from StandardAsyncExecutionActor
                             dockerStderr: String,                          // Calculated from StandardAsyncExecutionActor
                             inputs: Set[AwsBatchInput],
                             outputs: Set[AwsBatchFileOutput],
                             jobPaths: JobPaths,                    // Based on config, calculated in Job Paths, key to all things outside container
                             parameters: Seq[AwsBatchParameter]
                             ) {

  val Log = LoggerFactory.getLogger(AwsBatchJob.getClass)
  // TODO: Auth, endpoint
  lazy val client = BatchClient.builder()
                   // .credentialsProvider(...)
                   // .endpointOverride(...)
                   .build
  lazy val logsclient = CloudWatchLogsClient.builder()
                   // .credentialsProvider(...)
                   // .endpointOverride(...)
                   .build

  lazy val reconfiguredScript = {
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
    |Content-Type: multipart/alternative; boundary="${boundary}"

    |--${boundary}
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="rc.txt"
    |"
    |cat ${dockerRc}
    |echo "--${boundary}
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stdout.txt"
    |"
    |cat ${dockerStdout}
    |echo "--${boundary}
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stderr.txt"
    |"
    |cat ${dockerStderr}
    |echo "--${boundary}--"
    |exit $$(cat ${dockerRc})
    """).stripMargin
  }
  def submitJob(): Try[SubmitJobResponse] = Try {
    val taskId = jobDescriptor.key.call.fullyQualifiedName + "-" + jobDescriptor.key.index + "-" + jobDescriptor.key.attempt
    val workflow = jobDescriptor.workflowDescriptor
    val uniquePath = workflow.callable.name + "/" +
                     jobDescriptor.taskCall.callable.name + "/" +
                     workflow.id + "/" +
                     jobDescriptor.key.index + "/" +
                     jobDescriptor.key.attempt
    Log.info(s"""Submitting job to AWS Batch""")
    Log.info(s"""dockerImage: ${runtimeAttributes.dockerImage}""")
    Log.info(s"""jobQueueArn: ${runtimeAttributes.queueArn}""")
    Log.info(s"""taskId: $taskId""")
    Log.info(s"""hostpath root: $uniquePath""")

    // Build the Job definition before we submit. Eventually this should be
    // done separately and cached.
    val definitionArn = createDefinition(s"""${workflow.callable.name}-${jobDescriptor.taskCall.callable.name}""", uniquePath)

    var job : SubmitJobResponse = null // Use of null here due to underlying submitJob call
    for (inx <- 1 to 6) {
      try{
        if (job == null) {
          job = client.submitJob(SubmitJobRequest.builder()
                      .jobName(sanitize(jobDescriptor.taskCall.fullyQualifiedName))
                      .parameters(parameters.collect({ case i: AwsBatchInput => i.toStringString }).toMap.asJava)
                      .jobQueue(runtimeAttributes.queueArn)
                      .jobDefinition(definitionArn).build)
        }
      } catch {
        case e: ClientException if e.statusCode() == 404 && inx < 6 => {
          // RegisterJobDefinition is eventually consistent, so it may not be there
          Log.warn(s"Job Definition does not yet exist - retrying (${inx}/5)")
          Thread.sleep(100L * (inx ^ 2L))
        }
      }
    }
    job
  }

  /** Creates a job definition in AWS Batch
   *
   *  @param name Job definition name
   *  @return Arn for newly created job definition
   *
   */
  private def createDefinition(name: String, taskId: String): String = {
    val jobDefinitionBuilder = StandardAwsBatchJobDefinitionBuilder
    val jobDefinitionContext = AwsBatchJobDefinitionContext(runtimeAttributes,
                                                            taskId,
                                                            reconfiguredScript,
                                                            dockerRc,
                                                            dockerStdout,
                                                            dockerStderr,
                                                            jobDescriptor,
                                                            jobPaths,
                                                            inputs,
                                                            outputs)
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

    var arn = "" // Use of null here due to underlying submitJob call
    for (inx <- 1 to 6) {
      try {
        if (arn == "") {
          arn = client.registerJobDefinition(definitionRequest).jobDefinitionArn
        }
      } catch {
        case e: ClientException if e.statusCode() == 404 && inx < 6 => {
          // RegisterJobDefinition throws 404s every once in a while
          Log.warn(s"RegisterJobDefinition 404 - retrying (${inx}/5)")
          Thread.sleep(100L * (inx ^ 2L))
        }
        case e: ClientException if e.statusCode() == 409 => {
          // This could be a problem here, as we might have registerJobDefinition
          // but not describeJobDefinitions permissions. We've changed this to a
          // warning as a result of this potential
          Log.warn("Job definition already exists. Performing describe and retrieving latest revision")
          val defs = client
            .describeJobDefinitions(DescribeJobDefinitionsRequest.builder().jobDefinitionName(sanitize(name)).build())
            .jobDefinitions()
          // We need latest revision
          defs.get(defs.size() - 1).jobDefinitionArn()
        }
      }
    }
    arn
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
    name
      .replaceAll("[^A-Za-z0-9_\\-]", "_")
      .slice(0,128)

  /** Gets the status of a job by its Id, converted to a RunStatus
   *
   *  @param jobId Job ID as defined in AWS Batch
   *  @return Current RunStatus
   *
   */
  def status(jobId: String): Try[RunStatus] = {
     RunStatus.fromJobStatus(detail(jobId).status, jobId)
  }

  def detail(jobId: String): JobDetail = {
     val describeJobsResponse = client.describeJobs(DescribeJobsRequest.builder.jobs(jobId).build)
     val jobs = describeJobsResponse.jobs.asScala.toSeq
     jobs(0)
  }
  def rc(detail: JobDetail): Integer = {
     detail.container.exitCode
  }

  def output(detail: JobDetail): String = {
     // List<OutputEvent>
     val events = logsclient.getLogEvents(GetLogEventsRequest.builder
                                            // http://aws-java-sdk-javadoc.s3-website-us-west-2.amazonaws.com/latest/software/amazon/awssdk/services/batch/model/ContainerDetail.html#logStreamName--
                                            .logGroupName("/aws/batch/job")
                                            .logStreamName(detail.container.logStreamName)
                                            .startFromHead(true)
                                            .build).events.asScala.toSeq
     val eventMessages = for ( event <- events ) yield event.message
     eventMessages mkString "\n"
  }

  def abort(jobId: String): CancelJobResponse = {
    client.cancelJob(CancelJobRequest.builder.jobId(jobId).reason("cromwell abort called").build)
  }
}
