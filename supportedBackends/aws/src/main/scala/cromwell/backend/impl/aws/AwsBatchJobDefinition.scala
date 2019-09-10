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

import scala.language.postfixOps
import scala.collection.mutable.ListBuffer
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.JobPaths
import software.amazon.awssdk.services.batch.model.{ContainerProperties, KeyValuePair}
import wdl4s.parser.MemoryUnit
import cromwell.backend.impl.aws.io.AwsBatchVolume

import scala.collection.JavaConverters._
import java.io.ByteArrayOutputStream
import java.util.zip.GZIPOutputStream
import com.google.common.io.BaseEncoding
/*
 * This whole file might seem like a *lot* of ceremony just to build
 * the container properties. All this ceremony has been put in place for the
 * following two reasons:
 *
 * 1. It is consistent with other back ends
 * 2. It will allow for the possibility of builders other than the
 *    StandardAwsBatchJobDefinitionBuilder, which may prove useful in the
 *    future.
 *
 */
sealed trait AwsBatchJobDefinition {
  def containerProperties: ContainerProperties
}

trait AwsBatchJobDefinitionBuilder {

  /** Gets a builder, seeded with appropriate portions of the container properties
   *
   *  @param commandLine command line to execute within the container. Will be run in a shell context
   *  @param dockerImage docker image with which to run
   *  @return ContainerProperties builder ready for modification
   *
   */
  def builder(dockerImage: String): ContainerProperties.Builder =
    ContainerProperties.builder().image(dockerImage)

  def buildKVPair(key: String, value: String): KeyValuePair =
    KeyValuePair.builder.name(key).value(value).build


  def buildResources(builder: ContainerProperties.Builder,
                     context: AwsBatchJobDefinitionContext): ContainerProperties.Builder = {
    // The initial buffer should only contain one item - the hostpath of the
    // local disk mount point, which will be needed by the docker container
    // that copies data around
    val outputinfo = context.outputs.map(o => "%s,%s,%s,%s".format(o.name, o.s3key, o.local, o.mount))
      .mkString(";")
    val inputinfo = context.inputs.collect{case i: AwsBatchFileInput => i}
      .map(i => "%s,%s,%s,%s".format(i.name, i.s3key, i.local, i.mount))
      .mkString(";")

    val environment =
      context.runtimeAttributes.disks.collect{
        case d if d.name == "local-disk" =>  // this has s3 fiel system, needs all the env for the ecs-proxy
          List(buildKVPair("AWS_CROMWELL_LOCAL_DISK", d.mountPoint.toString),
          buildKVPair("AWS_CROMWELL_PATH",context.uniquePath),
          buildKVPair("AWS_CROMWELL_RC_FILE",context.dockerRcPath),
          buildKVPair("AWS_CROMWELL_STDOUT_FILE",context.dockerStdoutPath),
          buildKVPair("AWS_CROMWELL_STDERR_FILE",context.dockerStderrPath),
          buildKVPair("AWS_CROMWELL_CALL_ROOT",context.jobPaths.callExecutionRoot.toString),
          buildKVPair("AWS_CROMWELL_WORKFLOW_ROOT",context.jobPaths.workflowPaths.workflowRoot.toString),
          gzipKeyValuePair("AWS_CROMWELL_INPUTS", inputinfo),
          buildKVPair("AWS_CROMWELL_OUTPUTS",outputinfo))
      }.flatten

    def getVolPath(d:AwsBatchVolume) : Option[String] =  {
        d.fsType match {
          case "efs"  => None
          case _ =>  Option(context.uniquePath)
        }
    }

    builder
      .command(packCommand("/bin/bash", "-c", context.commandText).asJava)
      .memory(context.runtimeAttributes.memory.to(MemoryUnit.MB).amount.toInt)
      .vcpus(context.runtimeAttributes.cpu##)
      .volumes(context.runtimeAttributes.disks.map(d => d.toVolume(getVolPath(d))).asJava)
      .mountPoints(context.runtimeAttributes.disks.map(_.toMountPoint).asJava)
      .environment(environment.asJava)
  }

  private def gzip(data: String): String = {
    val byteArrayOutputStream = new ByteArrayOutputStream()
    val gzipOutputStream = new GZIPOutputStream(byteArrayOutputStream)
    gzipOutputStream.write(data.getBytes("UTF-8"))
    gzipOutputStream.close()

    BaseEncoding.base64().encode(byteArrayOutputStream.toByteArray())
  }

  private def packCommand(shell: String, options: String, mainCommand: String): Seq[String] = {
    val rc = new ListBuffer[String]()
    val lim = 1024
    val packedCommand = mainCommand.length() match {
      case len if len <= lim => mainCommand
      case len if len > lim => {
        rc += "gzipdata" // This is hard coded in our agent and must be the first item
        gzip(mainCommand)
      }
    }
    rc += shell
    rc += options
    rc += packedCommand
  }

  private def gzipKeyValuePair(prefix: String, data: String): KeyValuePair = {
    // This limit is somewhat arbitrary. The idea here is that if there are
    // a reasonable amount of inputs, we'll just stick with "AWS_CROMWELL_INPUTS".
    // This allows for easier user debugging as everything is in plain text.
    // However, if there are a ton of inputs, then a) the user doesn't want
    // to wade through all that muck when debugging, and b) we want to provide
    // as much room for the rest of the container overrides (of which there
    // are currently none, but this function doesn't know that).
    //
    // This whole thing is setup because there is an undisclosed AWS Batch
    // limit of 8k for container override JSON when it hits the endpoint.
    // In the few circumstances where there are a ton of inputs, we can
    // hit this limit. This particular value is highly compressible, though,
    // so we get a ton of runway this way. Incidentally, the overhead of
    // gzip+base64 can be larger on very small inputs than the input data
    // itself, so we don't want it super-small either. The proxy container
    // is designed to handle either of these variables.
    val lim = 512 // This is completely arbitrary and I think 512 may even be
                  // too big when looking at the value manually in the console
    data.length() match {
      case len if len <= lim => buildKVPair(prefix, data)
      case len if len > lim => buildKVPair(prefix + "_GZ", gzip(data))
    }
  }

  def build(context: AwsBatchJobDefinitionContext): AwsBatchJobDefinition
}

object StandardAwsBatchJobDefinitionBuilder extends AwsBatchJobDefinitionBuilder {
  def build(context: AwsBatchJobDefinitionContext): AwsBatchJobDefinition = {
    val builderInst = builder(context.runtimeAttributes.dockerImage)
    buildResources(builderInst, context)
    new StandardAwsBatchJobDefinitionBuilder(builderInst.build)
  }
}

case class StandardAwsBatchJobDefinitionBuilder private(containerProperties: ContainerProperties) extends AwsBatchJobDefinition

object AwsBatchJobDefinitionContext

case class AwsBatchJobDefinitionContext(
            runtimeAttributes: AwsBatchRuntimeAttributes,
            uniquePath: String,
            commandText: String,
            dockerRcPath: String,
            dockerStdoutPath: String,
            dockerStderrPath: String,
            jobDescriptor: BackendJobDescriptor,
            jobPaths: JobPaths,
            inputs: Set[AwsBatchInput],
            outputs: Set[AwsBatchFileOutput])
