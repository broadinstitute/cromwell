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

import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

class AwsBatchJobSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {

  behavior of "AwsBatchJob"

  it should "alter the commandScript for mime output" in {
    val script = """
    |tmpDir=mkdir -p "/cromwell-aws/cromwell-execution/wf_hello/2422ea26-2578-48b0-86e9-50cbdda7d70a/call-hello/tmp.39397e83" && echo "/cromwell-aws/cromwell-execution/wf_hello/2422ea26-2578-48b0-86e9-50cbdda7d70a/call-hello/tmp.39397e83"
    |chmod 777 "$tmpDir"
    |export _JAVA_OPTIONS=-Djava.io.tmpdir="$tmpDir"
    |export TMPDIR="$tmpDir"
    |export HOME="$HOME"
    |(
    |  cd /cromwell_root

    |)
    |(
    |  cd /cromwell_root


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


    |)
    |mv /cromwell_root/hello-rc.txt.tmp /cromwell_root/hello-rc.txt"""

    val boundary = "d283898f11be0dee6f9ac01470450bee"
    val expectedscript = """#!/bin/bash
    |mkdir -p /cromwell_root
    """.stripMargin + script.concat(s"""
    |echo "MIME-Version: 1.0
    |Content-Type: multipart/alternative; boundary="${boundary}"

    |--${boundary}
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="rc.txt"
    |"
    |cat /cromwell_root/hello-rc.txt
    |echo "--${boundary}
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stdout.txt"
    |"
    |cat /cromwell_root/hello-stdout.log
    |echo "--${boundary}
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename="stderr.txt"
    |"
    |cat /cromwell_root/hello-stderr.log
    |echo "--${boundary}--"
    """).stripMargin
    val job = AwsBatchJob(null, null, "commandLine", script,
      "/cromwell_root/hello-rc.txt", "/cromwell_root/hello-stdout.log", "/cromwell_root/hello-stderr.log", Seq.empty[AwsBatchParameter])

    job.reconfiguredScript should be(expectedscript)
  }

  it should "properly parse MIME output" in {
    val rawText = """MIME-Version: 1.0
    |Content-Type: multipart/alternative; boundary=d283898f11be0dee6f9ac01470450bee
    |--d283898f11be0dee6f9ac01470450bee
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename=rc.txt
    |0
    |--d283898f11be0dee6f9ac01470450bee
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename=stdout.txt
    |--d283898f11be0dee6f9ac01470450bee
    |Content-Type: text/plain
    |Content-Disposition: attachment; filename=stderr.txt
    |Hello World! Welcome to Cromwell . . . on AWS!
    |--d283898f11be0dee6f9ac01470450bee--""".stripMargin
    val (rc, stdout, stderr) = AwsBatchJob.parseOutput(rawText)
    rc should be (0)
    stdout should be ("")
    stderr should be ("Hello World! Welcome to Cromwell . . . on AWS!")
  }
}
