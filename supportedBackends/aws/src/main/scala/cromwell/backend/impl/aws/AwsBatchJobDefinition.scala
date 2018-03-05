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
import software.amazon.awssdk.services.batch.model.{
                                        // KeyValuePair,
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

  def buildResources(builder: ContainerProperties.Builder, runtimeAttributes: AwsBatchRuntimeAttributes): ContainerProperties.Builder = {
    builder
      .memory(runtimeAttributes.memory.toMegabytes.toInt)
      .vcpus(runtimeAttributes.cpu##)
      .volumes(runtimeAttributes.disks.map(_.toVolume).asJava)
      .mountPoints(runtimeAttributes.disks.map(_.toMountPoint).asJava)
  }

  def build(commandLine: String, runtimeAttributes: AwsBatchRuntimeAttributes, docker: String): AwsBatchJobDefinition
}

object StandardAwsBatchJobDefinitionBuilder extends AwsBatchJobDefinitionBuilder {
  def build(commandLine: String, runtimeAttributes: AwsBatchRuntimeAttributes, dockerImage: String): AwsBatchJobDefinition = {
    val builderInst = builder(commandLine, dockerImage)
    buildResources(builderInst, runtimeAttributes)
    new StandardAwsBatchJobDefinitionBuilder(builderInst.build)
  }
}

case class StandardAwsBatchJobDefinitionBuilder private(containerProperties: ContainerProperties) extends AwsBatchJobDefinition
