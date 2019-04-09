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

import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import common.collections.EnhancedCollections._
import cromwell.backend.BackendSpec
import cromwell.backend.io.JobPathsSpecHelper._
// import cromwell.cloudsupport.gcp.auth.AwsBatchAuthModeSpec
import cromwell.core.Tags.AwsTest
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsString}

class AwsBatchCallPathsSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {

  import BackendSpec._
  import AwsBatchTestConfig._

  behavior of "AwsBatchCallPaths"

  it should "map the correct filenames" taggedAs AwsTest in {
    // AwsBatchAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val configuration = new AwsBatchConfiguration(AwsBatchBackendConfigurationDescriptor)
    val workflowPaths = AwsBatchWorkflowPaths(
      workflowDescriptor,
      AnonymousCredentialsProvider.create.resolveCredentials(),
      configuration
    )

    val callPaths = AwsBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.returnCodeFilename should be("hello-rc.txt")
    callPaths.stderr.getFileName.pathAsString should be("s3://my-cromwell-workflows-bucket/hello-stderr.log")
    callPaths.stdout.getFileName.pathAsString should be("s3://my-cromwell-workflows-bucket/hello-stdout.log")
  }

  it should "map the correct paths" taggedAs AwsTest in {
    // AwsBatchAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val configuration = new AwsBatchConfiguration(AwsBatchBackendConfigurationDescriptor)
    val workflowPaths = AwsBatchWorkflowPaths(
      workflowDescriptor,
      AnonymousCredentialsProvider.create.resolveCredentials(),
      configuration
    )

    val callPaths = AwsBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.returnCode.pathAsString should
      be(s"s3://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-rc.txt")
    callPaths.stdout.pathAsString should
      be(s"s3://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stdout.log")
    callPaths.stderr.pathAsString should
      be(s"s3://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stderr.log")
  }

  it should "map the correct call context" taggedAs AwsTest in {
    // AwsBatchAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val configuration = new AwsBatchConfiguration(AwsBatchBackendConfigurationDescriptor)
    val workflowPaths = AwsBatchWorkflowPaths(
      workflowDescriptor,
      AnonymousCredentialsProvider.create.resolveCredentials(),
      configuration
    )

    val callPaths = AwsBatchJobPaths(workflowPaths, jobDescriptorKey)
    
    callPaths.callContext.root.pathAsString should
      be(s"s3://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should
      be(s"s3://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stdout.log")
    callPaths.callContext.stderr should
      be(s"s3://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stderr.log")
  }

}
