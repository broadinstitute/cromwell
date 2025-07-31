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

import scala.collection.mutable.ListBuffer
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.JobPaths
import cromwell.core.WorkflowOptions
import software.amazon.awssdk.services.batch.model.{
  ContainerProperties,
  EvaluateOnExit,
  Host,
  KeyValuePair,
  LinuxParameters,
  LogConfiguration,
  MountPoint,
  ResourceRequirement,
  ResourceType,
  RetryAction,
  RetryStrategy,
  Ulimit,
  Volume
}
import cromwell.backend.impl.aws.io.AwsBatchVolume

import scala.jdk.CollectionConverters._
import java.security.MessageDigest
import org.apache.commons.lang3.builder.{ToStringBuilder, ToStringStyle}
import org.slf4j.{Logger, LoggerFactory}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import scala.util.Success

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
  def retryStrategy: RetryStrategy
  def name: String

  override def toString: String =
    new ToStringBuilder(this, ToStringStyle.JSON_STYLE)
      .append("name", name)
      .append("containerProperties", containerProperties)
      .append("retryStrategy", retryStrategy)
      .build
}

trait AwsBatchJobDefinitionBuilder {
  val Log: Logger = LoggerFactory.getLogger(StandardAwsBatchJobDefinitionBuilder.getClass)

  /** Gets a builder, seeded with appropriate portions of the container properties
   *
   *  @param context AwsBatchJobDefinitionContext with all the runtime attributes
   *  @return ContainerProperties builder ready for modification and name
   *
   */
  def containerPropertiesBuilder(context: AwsBatchJobDefinitionContext): (ContainerProperties.Builder, String) = {

    def buildVolumes(disks: Seq[AwsBatchVolume], fsx: Option[List[String]]): List[Volume] = {

      val fsx_volumes = fsx.isDefined match {
        case true =>
          fsx.get.map(mnt => Volume.builder().name(mnt).host(Host.builder().sourcePath(s"/$mnt").build()).build())
        case false => List()
      }

      // all the configured disks plus the fetch and run volume and the aws-cli volume
      disks.map(d => d.toVolume()).toList ++ List(
        Volume
          .builder()
          .name("fetchAndRunScript")
          .host(Host.builder().sourcePath("/usr/local/bin/fetch_and_run.sh").build())
          .build(),
        // the aws-cli location on the EC2
        Volume
          .builder()
          .name("awsCliHome")
          .host(Host.builder().sourcePath("/usr/local/aws-cli").build())
          .build(),
        // the location of the instance-id on the host (set by cloud-init)
        // https://cloudinit.readthedocs.io/en/23.2.2/reference/faq.html#data
        Volume
          .builder()
          .name("instanceId")
          .host(Host.builder().sourcePath("/var/lib/cloud/data/instance-id").build())
          .build()
      ) ++ fsx_volumes
    }

    def buildMountPoints(disks: Seq[AwsBatchVolume], fsx: Option[List[String]]): List[MountPoint] = {

      val fsx_disks = fsx.isDefined match {
        case true =>
          fsx.get.map(mnt => MountPoint.builder().readOnly(false).sourceVolume(mnt).containerPath(s"/$mnt").build())
        case false => List()
      }

      // all the configured disks plus the fetch and run mount point and the AWS cli mount point
      disks.map(_.toMountPoint).toList ++ List(
        MountPoint
          .builder()
          .readOnly(true)
          .sourceVolume("fetchAndRunScript")
          .containerPath("/var/scratch/fetch_and_run.sh")
          .build(),
        MountPoint
          .builder()
          .readOnly(true)
          .sourceVolume("awsCliHome")
          // where the aws-cli will be on the container
          .containerPath("/usr/local/aws-cli")
          .build(),
        // the location of the instance-id on the container, used to tag the instance
        MountPoint
          .builder()
          .readOnly(true)
          .sourceVolume("instanceId")
          .containerPath("/var/lib/cloud/data/instance-id")
          .build()
      ) ++ fsx_disks
    }

    def buildUlimits(ulimits: Seq[Map[String, String]]): List[Ulimit] =
      ulimits
        .filter(_.nonEmpty)
        .map(u =>
          Ulimit
            .builder()
            .name(u("name"))
            .softLimit(u("softLimit").toInt)
            .hardLimit(u("hardLimit").toInt)
            .build()
        )
        .toList

    def buildName(imageName: String,
                  packedCommand: String,
                  volumes: List[Volume],
                  mountPoints: List[MountPoint],
                  env: Seq[KeyValuePair],
                  ulimits: List[Ulimit],
                  efsDelocalize: Boolean,
                  efsMakeMD5: Boolean,
                  tagResources: Boolean,
                  logGroupName: String,
                  sharedMemorySize: MemorySize,
                  fuseMount: Boolean,
                  jobTimeout: Int,
                  workflowOptions: WorkflowOptions
    ): String = {
      val roleArnStr = workflowOptions.getOrElse(AwsBatchWorkflowOptionKeys.JobRoleArn, "")
      s"$imageName:$packedCommand:${volumes.map(_.toString).mkString(",")}:${mountPoints.map(_.toString).mkString(",")}:${env
          .map(_.toString)
          .mkString(",")}:${ulimits.map(_.toString).mkString(",")}:${efsDelocalize.toString}:${efsMakeMD5.toString}:${tagResources.toString}:$logGroupName:${sharedMemorySize
          .to(MemoryUnit.MB)
          .amount
          .toInt}:${fuseMount.toString}:${jobTimeout}:$roleArnStr"
    }

    val environment = List.empty[KeyValuePair]
    val cmdName = context.runtimeAttributes.fileSystem match {
      case AWSBatchStorageSystems.s3 => "/var/scratch/fetch_and_run.sh"
      case _ => context.commandText
    }
    val packedCommand = packCommand("/bin/bash", "-c", cmdName)
    val volumes = buildVolumes(context.runtimeAttributes.disks, context.fsxMntPoint)
    val mountPoints = buildMountPoints(context.runtimeAttributes.disks, context.fsxMntPoint)
    val logGroupName = context.runtimeAttributes.logGroupName
    val logConfiguration = LogConfiguration
      .builder()
      .logDriver("awslogs")
      .options(
        Map(
          "awslogs-group" -> logGroupName
        ).asJava
      )
      .build()
    val ulimits = buildUlimits(context.runtimeAttributes.ulimits)
    val efsDelocalize = context.runtimeAttributes.efsDelocalize
    val efsMakeMD5 = context.runtimeAttributes.efsMakeMD5
    val tagResources = context.runtimeAttributes.tagResources

    val containerPropsName = buildName(
      context.runtimeAttributes.dockerImage,
      packedCommand.mkString(","),
      volumes,
      mountPoints,
      environment,
      ulimits,
      efsDelocalize,
      efsMakeMD5,
      tagResources,
      logGroupName,
      context.runtimeAttributes.sharedMemorySize,
      context.runtimeAttributes.fuseMount,
      context.runtimeAttributes.jobTimeout,
      context.workflowOptions
    )

    val linuxParametersBuilder = LinuxParameters
      .builder()
      .sharedMemorySize(context.runtimeAttributes.sharedMemorySize.to(MemoryUnit.MB).amount.toInt)

    if (context.runtimeAttributes.fuseMount) {
      linuxParametersBuilder
        .devices(
          List(
            software.amazon.awssdk.services.batch.model.Device
              .builder()
              .hostPath("/dev/fuse")
              .containerPath("/dev/fuse")
              .build()
          ).asJava
        )
    }

    val linuxParameters = linuxParametersBuilder.build()
    // simple true / false for now, depending on a single attribute
    val privileged = context.runtimeAttributes.fuseMount

    val builderWithBasicProperties = ContainerProperties
      .builder()
      .image(context.runtimeAttributes.dockerImage)
      .command(packedCommand.asJava)
      .resourceRequirements(
        ResourceRequirement
          .builder()
          .`type`(ResourceType.VCPU)
          .value(context.runtimeAttributes.cpu.##.toString)
          .build(),
        ResourceRequirement
          .builder()
          .`type`(ResourceType.MEMORY)
          .value(context.runtimeAttributes.memory.to(MemoryUnit.MB).amount.toInt.toString)
          .build()
      )
      .logConfiguration(logConfiguration)
      .volumes(volumes.asJava)
      .mountPoints(mountPoints.asJava)
      .environment(environment.asJava)
      .ulimits(ulimits.asJava)
      .linuxParameters(linuxParameters)
      .privileged(privileged)

    // Add job role ARN if specified
    val finalBuilder = context.workflowOptions.get(AwsBatchWorkflowOptionKeys.JobRoleArn) match {
      case Success(roleArn) => builderWithBasicProperties.jobRoleArn(roleArn)
      case _ => builderWithBasicProperties
    }

    (finalBuilder, containerPropsName)
  }

  def retryStrategyBuilder(context: AwsBatchJobDefinitionContext): (RetryStrategy.Builder, String) = {
    // We can add here the 'evaluateOnExit' statement
    var builder = RetryStrategy
      .builder()
      .attempts(context.runtimeAttributes.awsBatchRetryAttempts)

    var evaluations: Seq[EvaluateOnExit] = Seq()
    context.runtimeAttributes.awsBatchEvaluateOnExit.foreach { evaluate =>
      val evaluateBuilder = evaluate.foldLeft(EvaluateOnExit.builder()) { case (acc, (k, v)) =>
        (k.toLowerCase, v.toLowerCase) match {
          case ("action", "retry") => acc.action(RetryAction.RETRY)
          case ("action", "exit") => acc.action(RetryAction.EXIT)
          case ("onexitcode", _) => acc.onExitCode(v)
          case ("onreason", _) => acc.onReason(v)
          case ("onstatusreason", _) => acc.onStatusReason(v)
          case _ => acc
        }
      }
      evaluations = evaluations :+ evaluateBuilder.build()
    }

    builder = builder.evaluateOnExit(evaluations.asJava)

    (builder,
     s"${context.runtimeAttributes.awsBatchRetryAttempts.toString}${context.runtimeAttributes.awsBatchEvaluateOnExit.toString}"
    )
  }

  private def packCommand(shell: String, options: String, mainCommand: String): Seq[String] = {
    val rc = new ListBuffer[String]()
    val lim = 20480
    val packedCommand = mainCommand.length() match {
      case len if len <= lim => mainCommand
      case len if len > lim =>
        rc += "gzipdata" // This is hard coded in our agent and must be the first item
        gzip(mainCommand)
    }
    rc += shell
    rc += options
    rc += packedCommand
    rc.toList
  }

  def build(context: AwsBatchJobDefinitionContext): AwsBatchJobDefinition
}

object StandardAwsBatchJobDefinitionBuilder extends AwsBatchJobDefinitionBuilder {
  def build(context: AwsBatchJobDefinitionContext): AwsBatchJobDefinition = {

    val (containerPropsInst, containerPropsName) = containerPropertiesBuilder(context)
    val (retryStrategyInst, retryStrategyName) = retryStrategyBuilder(context)

    val name = buildName(context.runtimeAttributes.dockerImage, containerPropsName, retryStrategyName)

    new StandardAwsBatchJobDefinitionBuilder(containerPropsInst.build, retryStrategyInst.build, name)
  }

  def buildName(imageName: String, containerPropsName: String, retryStrategyName: String): String = {
    val str = s"$imageName:$containerPropsName:$retryStrategyName"

    val sha1 = MessageDigest
      .getInstance("SHA-1")
      .digest(str.getBytes("UTF-8"))
      .map("%02x".format(_))
      .mkString

    val prefix = s"cromwell_${imageName}_".slice(0, 88) // will be joined to a 40 character SHA1 for total length of 128

    sanitize(prefix + sha1)
  }
}

case class StandardAwsBatchJobDefinitionBuilder private (containerProperties: ContainerProperties,
                                                         retryStrategy: RetryStrategy,
                                                         name: String
) extends AwsBatchJobDefinition

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
  outputs: Set[AwsBatchFileOutput],
  fsxMntPoint: Option[List[String]],
  efsMntPoint: Option[String],
  efsMakeMD5: Option[Boolean],
  efsDelocalize: Option[Boolean],
  tagResources: Option[Boolean],
  workflowOptions: WorkflowOptions
) {

  override def toString: String =
    new ToStringBuilder(this, ToStringStyle.JSON_STYLE)
      .append("runtimeAttributes", runtimeAttributes)
      .append("commandText", commandText)
      .append("dockerRcPath", dockerRcPath)
      .append("dockerStderrPath", dockerStderrPath)
      .append("dockerStdoutPath", dockerStdoutPath)
      .append("jobDescriptor", jobDescriptor)
      .append("jobPaths", jobPaths)
      .append("inputs", inputs)
      .append("outputs", outputs)
      .append("fsxMntPoint", fsxMntPoint)
      .append("efsMntPoint", efsMntPoint)
      .append("efsMakeMD5", efsMakeMD5)
      .append("efsDelocalize", efsDelocalize)
      .append("tagResources", tagResources)
      .append("workflowOptions", workflowOptions)
      .build
}
