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
import cromwell.backend.BackendSpec._
import cromwell.backend.impl.aws.io.{AwsBatchJobPaths, AwsBatchWorkflowPaths, AwsBatchWorkingDisk}
import cromwell.backend.validation.ContinueOnReturnCodeFlag
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.TestKitSuite
import cromwell.core.path.DefaultPathBuilder
import cromwell.util.SampleWdl
import eu.timepit.refined.api.Refined
import eu.timepit.refined.auto._
import eu.timepit.refined.numeric._
import org.scalatest.PrivateMethodTester
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import software.amazon.awssdk.services.batch.model._
import spray.json.{JsObject, JsString}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.graph.CommandCallNode

import scala.jdk.javaapi.CollectionConverters

class AwsBatchJobSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with PrivateMethodTester {

  import AwsBatchTestConfig._

  System.setProperty("aws.region", "us-east-1")

  val script: String =
    """
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
  val jobDescriptor: BackendJobDescriptor =
    BackendJobDescriptor(workFlowDescriptor, jobKey, null, Map.empty, null, null, null)

  val jobPaths: AwsBatchJobPaths = AwsBatchJobPaths(workflowPaths, jobKey)
  val s3Inputs: Set[AwsBatchInput] = Set(
    AwsBatchFileInput("foo", "s3://bucket/foo", DefaultPathBuilder.get("foo"), AwsBatchWorkingDisk(), false, false)
  )
  val s3Outputs: Set[AwsBatchFileOutput] = Set(
    AwsBatchFileOutput("baa", "s3://bucket/somewhere/baa", DefaultPathBuilder.get("baa"), AwsBatchWorkingDisk(), false)
  )

  val cpu: Int Refined Positive = 2
  val sharedMemorySize: MemorySize = MemorySize(64, MemoryUnit.MB)
  val jobTimeout: Int = 3600

  val runtimeAttributes: AwsBatchRuntimeAttributes = new AwsBatchRuntimeAttributes(
    cpu = cpu,
    gpuCount = 0,
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
    awsBatchEvaluateOnExit = Vector(Map.empty[String, String]),
    ulimits = Vector(Map.empty[String, String]),
    efsDelocalize = false,
    efsMakeMD5 = false,
    fuseMount = false,
    fileSystem = "s3",
    sharedMemorySize = sharedMemorySize,
    jobTimeout = jobTimeout,
    logGroupName = "/aws/batch/job",
    additionalTags = Map("tag" -> "value")
  )

  val batchJobDefintion = AwsBatchJobDefinitionContext(
    runtimeAttributes = runtimeAttributes,
    commandText = "",
    dockerRcPath = "",
    dockerStdoutPath = "",
    dockerStderrPath = "",
    jobDescriptor = jobDescriptor,
    jobPaths = jobPaths,
    inputs = Set(),
    outputs = Set(),
    fsxMntPoint = None,
    None,
    None,
    None,
    None
  )

  val containerDetail: ContainerDetail = ContainerDetail.builder().exitCode(0).build()
  val jobDetail: JobDetail = JobDetail.builder().container(containerDetail).build

  private def generateBasicJob: AwsBatchJob = {
    val job = AwsBatchJob(
      jobDescriptor,
      runtimeAttributes,
      "commandLine",
      script,
      "/cromwell_root/hello-rc.txt",
      "/cromwell_root/hello-stdout.log",
      "/cromwell_root/hello-stderr.log",
      Seq.empty[AwsBatchInput].toSet,
      Seq.empty[AwsBatchFileOutput].toSet,
      jobPaths,
      Seq.empty[AwsBatchParameter],
      None,
      None,
      None,
      None,
      None,
      None,
      None,
      "",
      Map.empty
    )
    job
  }

  private def generateBasicJobForLocalFS: AwsBatchJob = {
    val job = AwsBatchJob(
      jobDescriptor,
      runtimeAttributes.copy(fileSystem = "local"),
      "commandLine",
      script,
      "/cromwell_root/hello-rc.txt",
      "/cromwell_root/hello-stdout.log",
      "/cromwell_root/hello-stderr.log",
      Seq.empty[AwsBatchInput].toSet,
      Seq.empty[AwsBatchFileOutput].toSet,
      jobPaths,
      Seq.empty[AwsBatchParameter],
      None,
      None,
      None,
      None,
      None,
      None,
      None,
      "",
      Map.empty
    )
    job
  }

  private def generateJobWithS3InOut: AwsBatchJob = {
    val job = AwsBatchJob(
      jobDescriptor,
      runtimeAttributes,
      "commandLine",
      script,
      "/cromwell_root/hello-rc.txt",
      "/cromwell_root/hello-stdout.log",
      "/cromwell_root/hello-stderr.log",
      s3Inputs,
      s3Outputs,
      jobPaths,
      Seq.empty[AwsBatchParameter],
      None,
      None,
      None,
      None,
      None,
      None,
      None,
      "",
      Map.empty
    )
    job
  }

  // TESTS BEGIN HERE
  behavior of "AwsBatchJob"

  it should "generate appropriate KV pairs for the container environment for S3" in {
    val job = generateBasicJob
    val generateEnvironmentKVPairs = PrivateMethod[List[KeyValuePair]](Symbol("generateEnvironmentKVPairs"))

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val kvPairs = job invokePrivate generateEnvironmentKVPairs("script-bucket", "prefix-", "key")
    kvPairs should contain(buildKVPair("BATCH_FILE_TYPE", "script"))
    kvPairs should contain(buildKVPair("BATCH_FILE_S3_URL", "s3://script-bucket/prefix-key"))
  }

  it should "generate appropriate KV pairs for the container environment for Local FS" in {
    val job = generateBasicJobForLocalFS
    val generateEnvironmentKVPairs = PrivateMethod[List[KeyValuePair]](Symbol("generateEnvironmentKVPairs"))

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val kvPairs = job invokePrivate generateEnvironmentKVPairs("script-bucket", "prefix-", "key")
    kvPairs should contain(buildKVPair("BATCH_FILE_TYPE", "script"))
    kvPairs should contain(buildKVPair("BATCH_FILE_S3_URL", ""))
  }

  it should "contain expected command script in reconfigured script" in {
    val job = generateBasicJob
    job.reconfiguredScript should include(script.replace("/cromwell_root", "/tmp/scratch"))
  }

  it should "add metadata environment variables to reconfigured script" in {
    val job = generateJobWithS3InOut
    job.reconfiguredScript should include("export AWS_METADATA_SERVICE_TIMEOUT=10\n")
    job.reconfiguredScript should include("export AWS_METADATA_SERVICE_NUM_ATTEMPTS=10\n")
  }

  it should "add s3 localize with retry function to reconfigured script" in {
    val job = generateBasicJob
    val retryFunctionText =
      s"""
         |function _s3_localize_with_retry() {
         |  local s3_path="$$1"
         |  # destination must be the path to a file and not just the directory you want the file in
         |  local destination="$$2"
         |  # if third option is specified, it is the optional tag (true / false)
         |  local is_optional="$${3:-false}"
         |  # if fourth option is specified, it is the locOptional tag (true / false)
         |  local loc_optional="$${4:-false}"
         |
         |  for i in {1..6};
         |  do
         |    # abort if tries are exhausted
         |    if [ "$$i" -eq 6 ]; then
         |        echo "failed to copy $$s3_path after $$(( $$i - 1 )) attempts."
         |        LOCALIZATION_FAILED=1
         |        return
         |    fi
         |    # check validity of source path
         |    if ! [[ "$$s3_path" =~ s3://([^/]+)/(.+) ]]; then
         |      echo "$$s3_path is not an S3 path with a bucket and key."
         |      LOCALIZATION_FAILED=1
         |      return
         |    fi
         |    ## if missing on s3 : check if optional:
         |    if ! /usr/local/aws-cli/v2/current/bin/aws s3 ls "$$s3_path" > /dev/null 2>&1 ; then
         |      if [[ "$$is_optional" == "true" ]]; then
         |        echo "Optional file '$$s3_path' does not exist. skipping localization"
         |      else
         |        echo "$$s3_path does not exist. skipping localization"
         |        LOCALIZATION_FAILED=1
         |      fi
         |      return
         |    fi
         |    # if localization is optional : skip
         |    if [[ "$$loc_optional" == "true" ]]; then
         |       echo "File $$s3_path does not have to be localized. Skipping localization"
         |       return 
         |    fi
         |    # copy
         |    /usr/local/aws-cli/v2/current/bin/aws s3 cp --no-progress "$$s3_path" "$$destination"  ||
         |        { echo "attempt $$i to copy $$s3_path failed" && sleep $$((7 * "$$i")) && continue; }
         |    # check data integrity
         |    _check_data_integrity "$$destination" "$$s3_path" ||
         |       { echo "data content length difference detected in attempt $$i to copy $$local_path failed" && sleep $$((7 * "$$i")) && continue; }
         |    # copy succeeded
         |    return
         |  done
         |}""".stripMargin

    job.reconfiguredScript should include(retryFunctionText)
  }

  it should "s3 delocalization with retry function in reconfigured script" in {
    val job = generateBasicJob
    val delocalizeText =
      s"""
         |function _s3_delocalize_with_retry() {
         |  # input variables
         |  local local_path="$$1"
         |  # destination must be the path to a file and not just the directory you want the file in
         |  local destination="$$2"
         |  # if third options is specified, it is the optional tag (true / false)
         |  local is_optional="$${3:-false}"
         |
         |  # get the multipart chunk size
         |  chunk_size=$$(_get_multipart_chunk_size "$$local_path")
         |  local MP_THRESHOLD=5368709120
         |  # then set them
         |  /usr/local/aws-cli/v2/current/bin/aws configure set default.s3.multipart_threshold $$MP_THRESHOLD
         |  /usr/local/aws-cli/v2/current/bin/aws configure set default.s3.multipart_chunksize $$chunk_size
         |
         |  # try & validate upload 5 times
         |  for i in {1..6};
         |  do
         |    # if tries exceeded : abort
         |    if [ "$$i" -eq 6 ]; then
         |        echo "failed to delocalize $$local_path after $$(( $$i - 1 )) attempts."
         |        DELOCALIZATION_FAILED=1
         |        return
         |    fi
         |    # if destination is not a bucket : abort
         |    if ! [[ "$$destination" =~ s3://([^/]+)/(.+) ]]; then
         |     echo "$$destination is not an S3 path with a bucket and key."
         |      DELOCALIZATION_FAILED=1
         |      return
         |    fi
         |    # copy ok or try again.
         |    if [[ -d "$$local_path" ]]; then
         |       # make sure to strip the trailing / in destination
         |       destination=$${destination%/}
         |       # glob directory. do recursive copy
         |       /usr/local/aws-cli/v2/current/bin/aws s3 cp --no-progress "$$local_path" "$$destination" --recursive --exclude "cromwell_glob_control_file" ||
         |         { echo "attempt $$i to copy globDir $$local_path failed" && sleep $$((7 * "$$i")) && continue; }
         |       # check integrity for each of the files (allow spaces)
         |       SAVEIFS="$$IFS"
         |       IFS=$$'\n'
         |       for FILE in $$(cd "$$local_path" ; ls | grep -v cromwell_glob_control_file); do
         |           _check_data_integrity "$$local_path/$$FILE" "$$destination/$$FILE" ||
         |               { echo "data content length difference detected in attempt $$i to copy $$local_path/$$FILE failed" && sleep $$((7 * "$$i")) && continue 2; }
         |       done
         |       IFS="$$SAVEIFS"
         |    # files : if exists or non-optional : must succeed
         |    elif [[ "$$is_optional" == "false" || -e "$$local_path" ]]; then
         |      /usr/local/aws-cli/v2/current/bin/aws s3 cp --no-progress "$$local_path" "$$destination" ||
         |         { echo "attempt $$i to copy $$local_path failed" && sleep $$((7 * "$$i")) && continue; }
         |      # check content length for data integrity
         |      _check_data_integrity "$$local_path" "$$destination" ||
         |         { echo "data content length difference detected in attempt $$i to copy $$local_path failed" && sleep $$((7 * "$$i")) && continue; }
         |    elif [[ "$$is_optional" == "true" && ! -e "$$local_path" ]]; then
         |      echo "Optional file '$$local_path' does not exist. skipping delocalization"
         |    # not optional, but missing : fail
         |    elif [[ "$$is_optional" == "false" && ! -e "$$local_path" ]]; then
         |      echo "$$local_path does not exist. skipping delocalization"
         |      DELOCALIZATION_FAILED=1
         |    fi
         |    # copy succeeded or not retrying
         |    return
         |  done
         |}""".stripMargin
    job.reconfiguredScript should include(delocalizeText)
  }

  it should "generate check data integrity in reconfigured script" in {
    val job = generateBasicJob
    val checkDataIntegrityBlock =
      s"""
         |function _check_data_integrity() {
         |  local local_path="$$1"
         |  local s3_path="$$2"
         |
         |  # remote : use content_length
         |  if [[ "$$s3_path" =~ s3://([^/]+)/(.+) ]]; then
         |        bucket="$${BASH_REMATCH[1]}"
         |        key="$${BASH_REMATCH[2]}"
         |  else
         |      # this is already checked in the caller function
         |      echo "$$s3_path is not an S3 path with a bucket and key."
         |      return 1
         |  fi
         |  s3_content_length=$$(/usr/local/aws-cli/v2/current/bin/aws s3api head-object --bucket "$$bucket" --key "$$key" --query 'ContentLength') ||
         |        { echo "Attempt to get head of object failed for $$s3_path." && return 1; }
         |  # local
         |  local_content_length=$$(LC_ALL=C ls -dnL -- "$$local_path" | awk '{print $$5; exit}' ) ||
         |        { echo "Attempt to get local content length failed for $$_local_path." && return 1; }
         |  # compare
         |  if [[ "$$s3_content_length" -eq "$$local_content_length" ]]; then
         |       true
         |  else
         |       false
         |  fi
         |}""".stripMargin
    job.reconfiguredScript should include(checkDataIntegrityBlock)
  }

  it should "generate get multipart chunk size in script" in {
    val job = generateBasicJob
    val getMultiplePartChunkSize =
      s"""
         |function _get_multipart_chunk_size() {
         |  local file_path="$$1"
         |  # missing files : skip.
         |  if [[ ! -e "$$file_path" ]]; then
         |    echo $$(( 5 * 1024 * 1024 ))
         |    return
         |  fi
         |  # file size
         |  file_size=$$(stat --printf="%s" "$$file_path")
         |  # chunk_size : you can have at most 10K parts with at least one 5MB part
         |  # this reflects the formula in s3-copy commands of cromwell (S3FileSystemProvider.java)
         |  #   => long partSize = Math.max((objectSize / 10000L) + 1, 5 * 1024 * 1024);
         |  a=$$(( ( file_size / 10000) + 1 ))
         |  b=$$(( 5 * 1024 * 1024 ))
         |  chunk_size=$$(( a > b ? a : b ))
         |  echo $$chunk_size
         |}
         |""".stripMargin

    job.reconfiguredScript should include(getMultiplePartChunkSize)
  }

  it should "generate postscript with output copy command in reconfigured script" in {
    val job = generateJobWithS3InOut
    val postscript =
      s"""
         |set -e
         |# (re-)add tags to include added volumes:
         |if [[ "false" == "true" ]]; then
         |  echo "*** TAGGING RESOURCES ***"
         |  _add_tags
         |fi
         |
         |echo '*** DELOCALIZING OUTPUTS ***'
         |DELOCALIZATION_FAILED=0
         |_s3_delocalize_with_retry "/tmp/scratch/baa" "s3://bucket/somewhere/baa" "false" 
         |
         |if [ -f "/tmp/scratch/hello-rc.txt" ]; then _s3_delocalize_with_retry "/tmp/scratch/hello-rc.txt" "${job.jobPaths.returnCode}" ; fi
         |if [ -f "/tmp/scratch/hello-stderr.log" ]; then _s3_delocalize_with_retry "/tmp/scratch/hello-stderr.log" "${job.jobPaths.standardPaths.error}"; fi
         |if [ -f "/tmp/scratch/hello-stdout.log" ]; then _s3_delocalize_with_retry "/tmp/scratch/hello-stdout.log" "${job.jobPaths.standardPaths.output}"; fi
         |
         |echo "DELOCALIZATION RESULT: $$DELOCALIZATION_FAILED"
         |if [[ $$DELOCALIZATION_FAILED -eq 1 ]]; then
         |  echo '*** DELOCALIZATION FAILED ***'
         |  echo '*** EXITING WITH RETURN CODE 1***'
         |  exit 1
         |else
         |  echo '*** COMPLETED DELOCALIZATION ***'
         |fi
         |echo '*** EXITING WITH RETURN CODE ***'
         |rc=$$(head -n 1 /tmp/scratch/hello-rc.txt)
         |echo $$rc
         |exit $$rc
         |}
         |""".stripMargin
    job.reconfiguredScript should include(postscript)
  }

  it should "generate preamble with input copy command in reconfigured script" in {
    val job = generateJobWithS3InOut
    val preamble =
      s"""
         |{
         |set -e
         |# tag instance and volumes to ensure tags are present in case of failure:
         |if [[ "false" == "true" ]]; then
         |  echo "*** TAGGING RESOURCES ***"
         |  _add_tags
         |fi
         |
         |echo '*** LOCALIZING INPUTS ***'
         |if [ ! -d /tmp/scratch ]; then mkdir /tmp/scratch && chmod 777 /tmp/scratch; fi
         |cd /tmp/scratch
         |# make sure localization completes successfully
         |LOCALIZATION_FAILED=0
         |_s3_localize_with_retry "s3://bucket/foo" "/tmp/scratch/foo" "false" "false" 
         |if [[ $$LOCALIZATION_FAILED -eq 1 ]]; then
         |  echo '*** LOCALIZATION FAILED ***'
         |  exit 1
         |else
         |  echo '*** COMPLETED LOCALIZATION ***'
         |fi
         |set +e
         |}
         |""".stripMargin

    job.reconfiguredScript should include(preamble)
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
    job.rc(jobDetail) should be(0)
  }

  it should "use RetryStrategy" in {
    val runtime = runtimeAttributes.copy(
      awsBatchEvaluateOnExit = Vector(Map("action" -> "EXIT", "onStatusReason" -> "Failed"))
    )

    val builder = RetryStrategy
      .builder()
      .attempts(1)
      .evaluateOnExit(
        EvaluateOnExit.builder().onStatusReason("Failed").action(RetryAction.EXIT).build()
      )
      .build()

    val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(batchJobDefintion.copy(runtimeAttributes = runtime))
    val expected = jobDefinition.retryStrategy
    expected should equal(builder)
  }

  it should "use RetryStrategy evaluateOnExit should be case insensitive" in {
    val runtime = runtimeAttributes.copy(
      awsBatchEvaluateOnExit = Vector(Map("aCtIoN" -> "EXIT", "onStatusReason" -> "Failed"))
    )

    val builder = RetryStrategy
      .builder()
      .attempts(1)
      .evaluateOnExit(
        EvaluateOnExit.builder().onStatusReason("Failed").action(RetryAction.EXIT).build()
      )
      .build()

    val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(batchJobDefintion.copy(runtimeAttributes = runtime))
    val expected = jobDefinition.retryStrategy
    expected should equal(builder)
  }

  it should "GPU is not set at job definition even if provided" in {
    val runtime = runtimeAttributes.copy(
      gpuCount = 1
    )

    val expected = List(
      ResourceRequirement.builder().`type`("VCPU").value(s"$cpu").build(),
      ResourceRequirement.builder().`type`("MEMORY").value("2048").build()
    )

    val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(batchJobDefintion.copy(runtimeAttributes = runtime))
    val actual = jobDefinition.containerProperties.resourceRequirements
    expected should equal(CollectionConverters.asScala(actual).toSeq)
  }

  it should "use default shared memory size of 64MB" in {
    val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(batchJobDefintion)
    val actual = jobDefinition.containerProperties.linuxParameters()
    val expected = LinuxParameters.builder().sharedMemorySize(64).build()
    expected should equal(actual)
  }

  it should "use user shared memory size if set" in {
    val runtime = runtimeAttributes.copy(
      gpuCount = 1,
      sharedMemorySize = MemorySize(100, MemoryUnit.MB)
    )
    val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(batchJobDefintion.copy(runtimeAttributes = runtime))
    val actual = jobDefinition.containerProperties.linuxParameters()
    val expected = LinuxParameters.builder().sharedMemorySize(100).build()
    expected should equal(actual)
  }
}

// ADD TEST FOR jobTimout
