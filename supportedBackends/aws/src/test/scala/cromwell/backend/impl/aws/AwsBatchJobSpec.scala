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

import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendJobDescriptorKey
import spray.json.{JsObject, JsString}
import common.collections.EnhancedCollections._

class AwsBatchJobSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {
  import AwsBatchTestConfig._

  behavior of "AwsBatchJob"

  /*
  DB 9/13/18
  This test is broken!  It does not pass *unless the user is logged in w/ aws cli* (i.e. "auth configure").  Our
  test server does not have such credentials and I was unsuccessful at setting this value manually.

  I'm merging this in the interest of time.

  It throws the following exception:
  Cause: java.lang.IllegalArgumentException: Could not build the path "s3://my-cromwell-workflows-bucket". It may refer to a filesystem not supported by this instance of Cromwell. Supported filesystems are: s3. Failures: s3: AWS region not provided (SdkClientException) Please refer to the documentation for more information on how to configure filesystems: http://cromwell.readthedocs.io/en/develop/backends/HPC/#filesystems
[info]   at cromwell.core.path.PathParsingException.<init>(PathParsingException.scala:5)
[info]   at cromwell.core.path.PathFactory$.$anonfun$buildPath$4(PathFactory.scala:64)
[info]   at scala.Option.getOrElse(Option.scala:121)
[info]   at cromwell.core.path.PathFactory$.buildPath(PathFactory.scala:58)
[info]   at cromwell.core.path.PathFactory.buildPath(PathFactory.scala:30)
[info]   at cromwell.core.path.PathFactory.buildPath$(PathFactory.scala:30)
[info]   at cromwell.backend.impl.aws.AwsBatchWorkflowPaths.buildPath(AwsBatchWorkflowPaths.scala:51)
[info]   at cromwell.backend.io.WorkflowPaths.executionRoot(WorkflowPaths.scala:34)
[info]   at cromwell.backend.io.WorkflowPaths.executionRoot$(WorkflowPaths.scala:34)
[info]   at cromwell.backend.impl.aws.AwsBatchWorkflowPaths.executionRoot$lzycompute(AwsBatchWorkflowPaths.scala:51)
   */
  it should "alter the commandScript for mime output" ignore {
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

    val boundary = "d283898f11be0dee6f9ac01470450bee"
    val expectedscript = script.concat(s"""
    |echo "MIME-Version: 1.0
    |Content-Type: multipart/alternative; boundary="$boundary"
    |
    |--$boundary
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="rc.txt"
    |"
    |cat /cromwell_root/hello-rc.txt
    |echo "--$boundary
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stdout.txt"
    |"
    |cat /cromwell_root/hello-stdout.log
    |echo "--$boundary
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stderr.txt"
    |"
    |cat /cromwell_root/hello-stderr.log
    |echo "--$boundary--"
    |exit $$(cat /cromwell_root/hello-rc.txt)
    """).stripMargin

    val workFlowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val configuration = new AwsBatchConfiguration(AwsBatchBackendConfigurationDescriptor)
    val workflowPaths = AwsBatchWorkflowPaths(
      workFlowDescriptor,
      AnonymousCredentialsProvider.create.resolveCredentials(),
      configuration
    )

    val call = workFlowDescriptor.callable.taskCallNodes.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val jobPaths = AwsBatchJobPaths(workflowPaths, jobKey)
    val job = AwsBatchJob(null, null, "commandLine", script,
      "/cromwell_root/hello-rc.txt", "/cromwell_root/hello-stdout.log", "/cromwell_root/hello-stderr.log", Seq.empty[AwsBatchInput].toSet, Seq.empty[AwsBatchFileOutput].toSet, jobPaths, Seq.empty[AwsBatchParameter], None)

    job.reconfiguredScript should be(expectedscript)
  }
}
