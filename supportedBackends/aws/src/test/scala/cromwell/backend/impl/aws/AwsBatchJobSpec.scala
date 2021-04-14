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

import common.collections.EnhancedCollections._
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.backend.BackendSpec._
import cromwell.backend.validation.ContinueOnReturnCodeFlag
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import eu.timepit.refined.api.Refined
import eu.timepit.refined.auto._
import eu.timepit.refined.numeric._
import org.scalatest.PrivateMethodTester
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import software.amazon.awssdk.services.batch.model.KeyValuePair
import spray.json.{JsObject, JsString}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.graph.CommandCallNode

class AwsBatchJobSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Mockito with PrivateMethodTester {
  import AwsBatchTestConfig._

  System.setProperty("aws.region", "us-east-1")

  val script = """
                 |tmpDir=mkdir -p "/cromwell-aws/cromwell-execution/wf_hello/2422ea26-2578-48b0-86e9-50cbdda7d70a/call-hello/tmp.39397e83" && echo "/cromwell-aws/cromwell-execution/wf_hello/2422ea26-2578-48b0-86e9-50cbdda7d70a/call-hello/tmp.39397e83"
                 |chmod 777 "$tmpDir"
                 |export _JAVA_OPTIONS=-Djava.io.tmpdir="$tmpDir"
                 |export TMPDIR="$tmpDir"
                 |export HOME="$HOME"
                 |(
                 |  cd /cromwell_root
                 |
                 |)
                 |(
                 |  cd /cromwell_root
                 |
                 |
                 |  echo "Hello World! Welcome to Cromwell . . . on AWS!" >&2
                 |)  > '/cromwell_root/hello-stdout.log' 2> '/cromwell_root/hello-stderr.log'
                 |echo $? > /cromwell_root/hello-rc.txt.tmp
                 |(
                 |  # add a .file in every empty directory to facilitate directory delocalization on the cloud
                 |  cd /cromwell_root
                 |  find . -type d -empty -print | xargs -I % touch %/.file
                 |)
                 |(
                 |  cd /cromwell_root
                 |  sync
                 |
                 |
                 |)
                 |mv /cromwell_root/hello-rc.txt.tmp /cromwell_root/hello-rc.txt"""

  val workFlowDescriptor: BackendWorkflowDescriptor = buildWdlWorkflowDescriptor(
    SampleWdl.HelloWorld.workflowSource(),
    inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
  )
  val configuration = new AwsBatchConfiguration(AwsBatchBackendConfigurationDescriptor)
  val workflowPaths: AwsBatchWorkflowPaths = AwsBatchWorkflowPaths(
    workFlowDescriptor,
    AnonymousCredentialsProvider.create,
    configuration
  )

  val call: CommandCallNode = workFlowDescriptor.callable.taskCallNodes.head
  val jobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)
  val jobPaths: AwsBatchJobPaths = AwsBatchJobPaths(workflowPaths, jobKey)

  val cpu: Int Refined Positive = 2
  val runtimeAttributes: AwsBatchRuntimeAttributes = new AwsBatchRuntimeAttributes(
      cpu = cpu,
      zones = Vector("us-east-1"),
      memory = MemorySize(2.0, MemoryUnit.GB),
      disks = Seq.empty,
      dockerImage = "ubuntu:latest",
      queueArn = "arn:aws:batch:us-east-1:123456789:job-queue/default-gwf-core",
      failOnStderr = true,
      continueOnReturnCode = ContinueOnReturnCodeFlag(false),
      noAddress = false,
      scriptS3BucketName = "script-bucket",
      awsBatchRetryAttempts = 1,
      ulimits = Vector(Map.empty[String, String]),
      fileSystem = "s3")

  private def generateBasicJob: AwsBatchJob = {
    val job = AwsBatchJob(null, runtimeAttributes, "commandLine", script,
      "/cromwell_root/hello-rc.txt", "/cromwell_root/hello-stdout.log", "/cromwell_root/hello-stderr.log",
      Seq.empty[AwsBatchInput].toSet, Seq.empty[AwsBatchFileOutput].toSet,
      jobPaths, Seq.empty[AwsBatchParameter], None)
    job
  }
  private def generateBasicJobForLocalFS: AwsBatchJob = {
    val job = AwsBatchJob(null, runtimeAttributes.copy(fileSystem="local"), "commandLine", script,
      "/cromwell_root/hello-rc.txt", "/cromwell_root/hello-stdout.log", "/cromwell_root/hello-stderr.log",
      Seq.empty[AwsBatchInput].toSet, Seq.empty[AwsBatchFileOutput].toSet,
      jobPaths, Seq.empty[AwsBatchParameter], None)
    job
  }

  // TESTS BEGIN HERE
  behavior of "AwsBatchJob"

  it should "have correctly named AWS constants" in {

    val job: AwsBatchJob = generateBasicJob

    job.AWS_RETRY_MODE should be ("AWS_RETRY_MODE")
    job.AWS_RETRY_MODE_DEFAULT_VALUE should be ("adaptive")
    job.AWS_MAX_ATTEMPTS should be ("AWS_MAX_ATTEMPTS")
    job.AWS_MAX_ATTEMPTS_DEFAULT_VALUE should be ("14")
  }

  it should "generate appropriate KV pairs for the container environment for S3" in {
    val job = generateBasicJob
    val generateEnvironmentKVPairs = PrivateMethod[List[KeyValuePair]]('generateEnvironmentKVPairs)

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val kvPairs = job invokePrivate generateEnvironmentKVPairs("script-bucket", "prefix-", "key")

    kvPairs should contain (buildKVPair(job.AWS_MAX_ATTEMPTS, job.AWS_MAX_ATTEMPTS_DEFAULT_VALUE))
    kvPairs should contain (buildKVPair(job.AWS_RETRY_MODE, "adaptive"))
    kvPairs should contain (buildKVPair("BATCH_FILE_TYPE", "script"))
    kvPairs should contain (buildKVPair("BATCH_FILE_S3_URL", "s3://script-bucket/prefix-key"))
  }

  it should "generate appropriate KV pairs for the container environment for Local FS" in {
    val job = generateBasicJobForLocalFS
    val generateEnvironmentKVPairs = PrivateMethod[List[KeyValuePair]]('generateEnvironmentKVPairs)

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val kvPairs = job invokePrivate generateEnvironmentKVPairs("script-bucket", "prefix-", "key")

    kvPairs should contain (buildKVPair(job.AWS_MAX_ATTEMPTS, job.AWS_MAX_ATTEMPTS_DEFAULT_VALUE))
    kvPairs should contain (buildKVPair(job.AWS_RETRY_MODE, "adaptive"))
    kvPairs should contain (buildKVPair("BATCH_FILE_TYPE", "script"))
    kvPairs should contain (buildKVPair("BATCH_FILE_S3_URL", ""))
  }
}
