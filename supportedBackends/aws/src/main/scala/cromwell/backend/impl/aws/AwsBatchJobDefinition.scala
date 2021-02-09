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
import software.amazon.awssdk.services.batch.model.{ContainerProperties, Host, KeyValuePair, MountPoint, Volume}
import cromwell.backend.impl.aws.io.AwsBatchVolume

import scala.collection.JavaConverters._
import java.security.MessageDigest

import org.apache.commons.lang3.builder.{ToStringBuilder, ToStringStyle}
import org.slf4j.{Logger, LoggerFactory}
import wdl4s.parser.MemoryUnit


/**
  * Responsible for the creation of the job definition.
  *
  * In AWS Batch, every job must have a definition.
  * In the current backend, the job definition is a 1:1 with a job.
  * This whole file might seem like a '''lot''' of ceremony just to build
  * the container properties. All this ceremony has been put in place for the
  * following two reasons:
  *
  *   1. It is consistent with other back ends
  *   2. It will allow for the possibility of builders other than the `StandardAwsBatchJobDefinitionBuilder`,
  *      which may prove useful in the future.
  *
  */
sealed trait AwsBatchJobDefinition {
  def containerProperties: ContainerProperties
  def name: String

  override def toString: String = {
    new ToStringBuilder(this, ToStringStyle.JSON_STYLE)
      .append("name", name)
      .append("containerProperties", containerProperties)
      .build
  }
}

trait AwsBatchJobDefinitionBuilder {
  val Log: Logger = LoggerFactory.getLogger(StandardAwsBatchJobDefinitionBuilder.getClass)


  /** Gets a builder, seeded with appropriate portions of the container properties
   *
   *  @param dockerImage docker image with which to run
   *  @return ContainerProperties builder ready for modification
   *
   */
  def builder(dockerImage: String): ContainerProperties.Builder =
    ContainerProperties.builder().image(dockerImage)


  def buildResources(builder: ContainerProperties.Builder,
                     context: AwsBatchJobDefinitionContext): (ContainerProperties.Builder, String) = {
    // The initial buffer should only contain one item - the hostpath of the
    // local disk mount point, which will be needed by the docker container
    // that copies data around

    val environment = List.empty[KeyValuePair]


    def buildVolumes(disks: Seq[AwsBatchVolume]): List[Volume] = {

      //all the configured disks plus the fetch and run volume and the aws-cli volume
      disks.map(d => d.toVolume()).toList ++ List(
        Volume.builder()
        .name("fetchAndRunScript")
        .host(Host.builder().sourcePath("/usr/local/bin/fetch_and_run.sh").build())
        .build(),
        //the aws-cli location on the EC2
        Volume.builder()
          .name("awsCliHome")
          .host(Host.builder().sourcePath("/usr/local/aws-cli").build())
          .build()
      )
    }

    def buildMountPoints(disks: Seq[AwsBatchVolume]): List[MountPoint] = {

      //all the configured disks plus the fetch and run mount point and the AWS cli mount point
      disks.map(_.toMountPoint).toList ++ List(
        MountPoint.builder()
          .readOnly(true)
          .sourceVolume("fetchAndRunScript")
          .containerPath("/var/scratch/fetch_and_run.sh")
          .build(),

        MountPoint.builder()
          .readOnly(true)
          .sourceVolume("awsCliHome")
          //where the aws-cli will be on the container
          .containerPath("/usr/local/aws-cli")
          .build()
      )
    }

    def buildName(imageName: String, packedCommand: String, volumes: List[Volume], mountPoints: List[MountPoint], env: Seq[KeyValuePair]): String = {
      val str = s"$imageName:$packedCommand:${volumes.map(_.toString).mkString(",")}:${mountPoints.map(_.toString).mkString(",")}:${env.map(_.toString).mkString(",")}"

      val sha1 = MessageDigest.getInstance("SHA-1")
            .digest( str.getBytes("UTF-8") )
            .map("%02x".format(_)).mkString

      val prefix = s"cromwell_$imageName".slice(0,88) // will be joined to a 40 character SHA1 for total length of 128

      sanitize(prefix + sha1)
    }


    val cmdName = context.runtimeAttributes.fileSystem match {
       case AWSBatchStorageSystems.s3 => "/var/scratch/fetch_and_run.sh"
       case _ =>  context.commandText
    }
    val packedCommand = packCommand("/bin/bash", "-c", cmdName)
    val volumes =  buildVolumes( context.runtimeAttributes.disks )
    val mountPoints = buildMountPoints( context.runtimeAttributes.disks)
    val jobDefinitionName = buildName(
      context.runtimeAttributes.dockerImage,
      packedCommand.mkString(","),
      volumes,
      mountPoints,
      environment
    )

    (builder
       .command(packedCommand.asJava)
        .memory(context.runtimeAttributes.memory.to(MemoryUnit.MB).amount.toInt)
        .vcpus(context.runtimeAttributes.cpu##)
        .volumes( volumes.asJava)
        .mountPoints( mountPoints.asJava)
        .environment(environment.asJava),

      jobDefinitionName)
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

  def build(context: AwsBatchJobDefinitionContext): AwsBatchJobDefinition
}

object StandardAwsBatchJobDefinitionBuilder extends AwsBatchJobDefinitionBuilder {
  def build(context: AwsBatchJobDefinitionContext): AwsBatchJobDefinition = {
    //instantiate a builder with the name of the docker image
    val builderInst = builder(context.runtimeAttributes.dockerImage)
    val (b, name) = buildResources(builderInst, context)

    new StandardAwsBatchJobDefinitionBuilder(b.build, name)
  }
}

case class StandardAwsBatchJobDefinitionBuilder private(containerProperties: ContainerProperties, name: String) extends AwsBatchJobDefinition

object AwsBatchJobDefinitionContext

case class AwsBatchJobDefinitionContext(
            runtimeAttributes: AwsBatchRuntimeAttributes,
            commandText: String,
            dockerRcPath: String,
            dockerStdoutPath: String,
            dockerStderrPath: String,
            jobDescriptor: BackendJobDescriptor,
            jobPaths: JobPaths,
            inputs: Set[AwsBatchInput],
            outputs: Set[AwsBatchFileOutput]){

  override def toString: String = {
    new ToStringBuilder(this, ToStringStyle.JSON_STYLE)
      .append("runtimeAttributes", runtimeAttributes)
      .append("commandText", commandText)
      .append("dockerRcPath", dockerRcPath)
      .append("dockerStderrPath",dockerStderrPath)
      .append("dockerStdoutPath", dockerStdoutPath)
      .append("jobDescriptor", jobDescriptor)
      .append("jobPaths", jobPaths)
      .append("inputs", inputs)
      .append("outputs", outputs)
      .build
  }
}
