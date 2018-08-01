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
import cromwell.backend.BackendJobDescriptor
import cromwell.core.path.Path
import software.amazon.awssdk.services.batch.model.{
                                        KeyValuePair,
                                        ContainerProperties //,
                                        // MountPoint,
                                        // Ulimit,
                                        // Volume,
                                      }

import scala.collection.JavaConverters._

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
  def builder(commandLine: String, dockerImage: String): ContainerProperties.Builder =
    ContainerProperties.builder().command("/bin/bash", "-c", commandLine).image(dockerImage)

  def buildKVPair(key: String, value: String): KeyValuePair =
    KeyValuePair.builder.name(key).value(value).build

  def buildResources(builder: ContainerProperties.Builder,
                     runtimeAttributes: AwsBatchRuntimeAttributes,
                     uniquePath: String,
                     rcPath: String,
                     jobDescriptor: BackendJobDescriptor,
                     callExecutionRoot: Path): ContainerProperties.Builder = {
    // The initial buffer should only contain one item - the hostpath of the
    // local disk mount point, which will be needed by the docker container
    // that copies data around
    val environment =
      runtimeAttributes.disks.collect{
        case d if d.name == "local-disk" =>
          buildKVPair("AWS_CROMWELL_LOCAL_DISK", d.mountPoint.toString)
      }.toBuffer
    environment.append(buildKVPair("AWS_CROMWELL_PATH",uniquePath))
    environment.append(buildKVPair("AWS_CROMWELL_RC_FILE",rcPath))
    environment.append(buildKVPair("AWS_CROMWELL_CALL_ROOT",callExecutionRoot.toString))
    // We want a marker for the process monitor to know for certain this container is from
    // Cromwell. Rather than relying on some of the other environment variables, we're using
    // a unique marker variable that doesn't mean much to anything else (and therefore has
    // no reason to be passed on. Also, by mispelling "CROMWELL" as "CRMWLL" unlike the other variables,
    // any process scanning for all "AWS_CROMWELL" will skip by this particular one.
    environment.append(buildKVPair("AWS_CRMWLL_PROCESS_MONITOR_MARKER","aws_batch_cromwell_process_monitor_marker"))

    builder
      .memory(runtimeAttributes.memory.toMegabytes.toInt)
      .vcpus(runtimeAttributes.cpu##)
      .volumes(runtimeAttributes.disks.map(_.toVolume(uniquePath)).asJava)
      .mountPoints(runtimeAttributes.disks.map(_.toMountPoint).asJava)
      .environment(environment.asJava)
  }

  def build(commandLine: String, runtimeAttributes: AwsBatchRuntimeAttributes,
            docker: String, uniquePath: String, rcPath: String, jobDescriptor: BackendJobDescriptor, callExecutionRoot: Path): AwsBatchJobDefinition
}

object StandardAwsBatchJobDefinitionBuilder extends AwsBatchJobDefinitionBuilder {
  def build(commandLine: String, runtimeAttributes: AwsBatchRuntimeAttributes,
            dockerImage: String, uniquePath: String, rcPath: String, jobDescriptor: BackendJobDescriptor, callExecutionRoot: Path): AwsBatchJobDefinition = {
    val builderInst = builder(commandLine, dockerImage)
    buildResources(builderInst, runtimeAttributes, uniquePath, rcPath, jobDescriptor, callExecutionRoot)
    new StandardAwsBatchJobDefinitionBuilder(builderInst.build)
  }
}

case class StandardAwsBatchJobDefinitionBuilder private(containerProperties: ContainerProperties) extends AwsBatchJobDefinition
