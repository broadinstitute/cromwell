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
import cromwell.backend.impl.aws.io.AwsBatchWorkingDisk
import cromwell.backend.validation.ContinueOnReturnCodeFlag
import cromwell.core.path.DefaultPathBuilder
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
import software.amazon.awssdk.services.batch.model.{ContainerDetail, JobDetail, KeyValuePair}
import spray.json.{JsObject, JsString}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.graph.CommandCallNode

class AwsBatchJobSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Mockito with PrivateMethodTester {
  import AwsBatchTestConfig._

  System.setProperty("aws.region", "us-east-1")

  val script: String = """
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
                 |mv /cromwell_root/hello-rc.txt.tmp /cromwell_root/hello-rc.txt""".stripMargin

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
  val s3Inputs: Set[AwsBatchInput] = Set(AwsBatchFileInput("foo", "s3://bucket/foo", DefaultPathBuilder.get("foo"), AwsBatchWorkingDisk()))
  val s3Outputs: Set[AwsBatchFileOutput] = Set(AwsBatchFileOutput("baa", "s3://bucket/somewhere/baa", DefaultPathBuilder.get("baa"), AwsBatchWorkingDisk()))

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

  val containerDetail: ContainerDetail = ContainerDetail.builder().exitCode(0).build()
  val jobDetail: JobDetail = JobDetail.builder().container(containerDetail).build

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
  private def generateJobWithS3InOut: AwsBatchJob = {
    val job = AwsBatchJob(null, runtimeAttributes, "commandLine", script,
      "/cromwell_root/hello-rc.txt", "/cromwell_root/hello-stdout.log", "/cromwell_root/hello-stderr.log",
      s3Inputs, s3Outputs,
      jobPaths, Seq.empty[AwsBatchParameter], None)
    job
  }

  // TESTS BEGIN HERE
  behavior of "AwsBatchJob"
  it should "generate appropriate KV pairs for the container environment for S3" in {
    val job = generateBasicJob
    val generateEnvironmentKVPairs = PrivateMethod[List[KeyValuePair]]('generateEnvironmentKVPairs)

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val kvPairs = job invokePrivate generateEnvironmentKVPairs("script-bucket", "prefix-", "key")
    kvPairs should contain (buildKVPair("BATCH_FILE_TYPE", "script"))
    kvPairs should contain (buildKVPair("BATCH_FILE_S3_URL", "s3://script-bucket/prefix-key"))
  }

  it should "generate appropriate KV pairs for the container environment for Local FS" in {
    val job = generateBasicJobForLocalFS
    val generateEnvironmentKVPairs = PrivateMethod[List[KeyValuePair]]('generateEnvironmentKVPairs)

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val kvPairs = job invokePrivate generateEnvironmentKVPairs("script-bucket", "prefix-", "key")
    kvPairs should contain (buildKVPair("BATCH_FILE_TYPE", "script"))
    kvPairs should contain (buildKVPair("BATCH_FILE_S3_URL", ""))
  }

  it should "contain expected command script in reconfigured script" in {
    val job = generateBasicJob
    job.reconfiguredScript should include (script.replace("/cromwell_root", "/tmp/scratch"))
  }

  it should "add metadata environment variables to reconfigured script" in {
    val job = generateJobWithS3InOut
    job.reconfiguredScript should include ("export AWS_METADATA_SERVICE_TIMEOUT=10\n")
    job.reconfiguredScript should include ("export AWS_METADATA_SERVICE_NUM_ATTEMPTS=10\n")
  }

  it should "add s3 localize with retry function to reconfigured script" in {
    val job = generateBasicJob
    val retryFunctionText = s"""
                              |function _s3_localize_with_retry() {
                              |  local s3_path=$$1
                              |  # destination must be the path to a file and not just the directory you want the file in
                              |  local destination=$$2
                              |
                              |  for i in {1..5};
                              |  do
                              |    if [[ $$s3_path =~ s3://([^/]+)/(.+) ]]; then
                              |        bucket="$${BASH_REMATCH[1]}"
                              |        key="$${BASH_REMATCH[2]}"
                              |        content_length=$$(/usr/local/aws-cli/v2/current/bin/aws  s3api head-object --bucket "$$bucket" --key "$$key" --query 'ContentLength')
                              |    else
                              |      echo "$$s3_path is not an S3 path with a bucket and key. aborting"
                              |      exit 1
                              |    fi
                              |    /usr/local/aws-cli/v2/current/bin/aws  s3 cp --no-progress "$$s3_path" "$$destination" &&
                              |    [[ $$(LC_ALL=C ls -dn -- "$$destination" | awk '{print $$5; exit}') -eq "$$content_length" ]] && break ||
                              |    echo "attempt $$i to copy $$s3_path failed";
                              |
                              |    if [ "$$i" -eq 5 ]; then
                              |        echo "failed to copy $$s3_path after $$i attempts. aborting"
                              |        exit 2
                              |    fi
                              |    sleep $$((7 * "$$i"))
                              |  done
                              |}
                              |""".stripMargin

    job.reconfiguredScript should include (retryFunctionText)
  }

  it should "generate postscript with output copy command in reconfigured script" in {
    val job = generateJobWithS3InOut
    val postscript =
      s"""
         |{
         |set -e
         |echo '*** DELOCALIZING OUTPUTS ***'
         |
         |/usr/local/aws-cli/v2/current/bin/aws  s3 cp --no-progress /tmp/scratch/baa s3://bucket/somewhere/baa
         |
         |
         |if [ -f /tmp/scratch/hello-rc.txt ]; then /usr/local/aws-cli/v2/current/bin/aws  s3 cp --no-progress /tmp/scratch/hello-rc.txt ${job.jobPaths.returnCode} ; fi
         |
         |if [ -f /tmp/scratch/hello-stderr.log ]; then /usr/local/aws-cli/v2/current/bin/aws  s3 cp --no-progress /tmp/scratch/hello-stderr.log ${job.jobPaths.standardPaths.error}; fi
         |if [ -f /tmp/scratch/hello-stdout.log ]; then /usr/local/aws-cli/v2/current/bin/aws  s3 cp --no-progress /tmp/scratch/hello-stdout.log ${job.jobPaths.standardPaths.output}; fi
         |
         |echo '*** COMPLETED DELOCALIZATION ***'
         |echo '*** EXITING WITH RETURN CODE ***'
         |rc=$$(head -n 1 /tmp/scratch/hello-rc.txt)
         |echo $$rc
         |exit $$rc
         |}
         |""".stripMargin
    job.reconfiguredScript should include (postscript)
  }

  it should "generate preamble with input copy command in reconfigured script" in {
    val job = generateJobWithS3InOut
    val preamble =
      s"""
         |{
         |set -e
         |echo '*** LOCALIZING INPUTS ***'
         |if [ ! -d /tmp/scratch ]; then mkdir /tmp/scratch && chmod 777 /tmp/scratch; fi
         |cd /tmp/scratch
         |_s3_localize_with_retry s3://bucket/foo /tmp/scratch/foo
         |echo '*** COMPLETED LOCALIZATION ***'
         |set +e
         |}
         |""".stripMargin

    job.reconfiguredScript should include (preamble)
  }

  it should "contain AWS Service clients" in {
    val job = generateBasicJob
    job.batchClient should not be null
    job.s3Client should not be null
    job.cloudWatchLogsClient should not be null
  }

  it should "have correct script prefix" in {
    val job = generateBasicJob
    job.scriptKeyPrefix should equal("scripts/")
  }

  it should "return correct RC code given Batch Job Detail" in {
    val containerDetail: ContainerDetail = ContainerDetail.builder().exitCode(0).build
    val jobDetail: JobDetail = JobDetail.builder().container(containerDetail).build
    val job = generateBasicJob
    job.rc(jobDetail) should be (0)
  }



}
