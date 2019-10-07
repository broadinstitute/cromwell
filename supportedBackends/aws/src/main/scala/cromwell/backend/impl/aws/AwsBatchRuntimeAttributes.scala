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

import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import cromwell.backend.impl.aws.io.{AwsBatchVolume, AwsBatchWorkingDisk}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{BooleanRuntimeAttributesValidation, _}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wom.types._
import wom.values._


case class AwsBatchRuntimeAttributes(cpu: Int Refined Positive,
                                zones: Vector[String],
                                memory: MemorySize,
                                disks: Seq[AwsBatchVolume],
                                dockerImage: String,
                                queueArn: String,
                                failOnStderr: Boolean,
                                continueOnReturnCode: ContinueOnReturnCode,
                                noAddress: Boolean,
                                     fileSystem:String= "s3")

object AwsBatchRuntimeAttributes {

  val QueueArnKey = "queueArn"

  val ZonesKey = "zones"
  private val ZonesDefaultValue = WomString("us-east-1a")

  val NoAddressKey = "noAddress"
  private val noAddressValidationInstance = new BooleanRuntimeAttributesValidation(NoAddressKey)
  private val NoAddressDefaultValue = WomBoolean(false)

  // TODO: Determine good volume format
  val DisksKey = "disks"
  private val DisksDefaultValue = WomString(s"${AwsBatchWorkingDisk.Name}")

  private val MemoryDefaultValue = "2 GB"

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation.instance
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def cpuMinValidation(runtimeConfig: Option[Config]):RuntimeAttributesValidation[Int Refined Positive] = CpuValidation.instanceMin
    .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def disksValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Seq[AwsBatchVolume]] = DisksValidation
    .withDefault(DisksValidation.configDefaultWomValue(runtimeConfig) getOrElse DisksDefaultValue)

  private def zonesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Vector[String]] = ZonesValidation
    .withDefault(ZonesValidation.configDefaultWomValue(runtimeConfig) getOrElse ZonesDefaultValue)

  private def memoryValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[MemorySize] = {
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.MemoryKey,
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, runtimeConfig) getOrElse MemoryDefaultValue)
  }

  private def memoryMinValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[MemorySize] = {
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.MemoryMinKey,
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryMinKey, runtimeConfig) getOrElse MemoryDefaultValue)
  }

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] = noAddressValidationInstance
    .withDefault(noAddressValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue)

  private val dockerValidation: RuntimeAttributesValidation[String] = DockerValidation.instance

  private def queueArnValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[String] =
    QueueArnValidation.withDefault(QueueArnValidation.configDefaultWomValue(runtimeConfig) getOrElse
      (throw new RuntimeException("queueArn is required")))

  def runtimeAttributesBuilder(configuration: AwsBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = configuration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      cpuValidation(runtimeConfig),
      cpuMinValidation(runtimeConfig),
      disksValidation(runtimeConfig),
      zonesValidation(runtimeConfig),
      memoryValidation(runtimeConfig),
      memoryMinValidation(runtimeConfig),
      noAddressValidation(runtimeConfig),
      dockerValidation,
      queueArnValidation(runtimeConfig)
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): AwsBatchRuntimeAttributes = {
    val cpu: Int Refined Positive = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val memory: MemorySize = RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val disks: Seq[AwsBatchVolume] = RuntimeAttributesValidation.extract(disksValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val docker: String = RuntimeAttributesValidation.extract(dockerValidation, validatedRuntimeAttributes)
    val queueArn: String = RuntimeAttributesValidation.extract(queueArnValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val failOnStderr: Boolean = RuntimeAttributesValidation.extract(failOnStderrValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val noAddress: Boolean = RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    new AwsBatchRuntimeAttributes(
      cpu,
      zones,
      memory,
      disks,
      docker,
      queueArn,
      failOnStderr,
      continueOnReturnCode,
      noAddress
    )
  }
}

object QueueArnValidation extends ArnValidation(AwsBatchRuntimeAttributes.QueueArnKey) {
  // queue arn format can be found here
  // https://docs.aws.amazon.com/en_us/general/latest/gr/aws-arns-and-namespaces.html#arn-syntax-batch
  // arn:aws:batch:region:account-id:job-queue/queue-name
  override protected val arnRegex =
  s"""
      (?x)                            # Turn on comments and whitespace insensitivity
      (arn)                           # Every AWS ARN starts with "arn"
      :
      (                               # Begin capturing ARN for partition
        aws                           # Required part of a partition
        (-[a-z]+){0,2}                # Optional part like "-cn" or "-us-gov". Therefore, the whole partition may look like "aws", "aws-cn", "aws-us-gov"
      )                               # End capturing ARN for partition
      :
      (batch)                         # Job queues is contained in AWS Batch
      :
      (                               # Begin capturing ARN for region
        [a-z]{2}                      # ALPHA-2 county code
        (-[a-z]+){1,2}                # Specific region like "-west" or "-southeast". In case of AWS GovCloud it may be like "-gov-east"
        -\\d                          # Sequence number of the region
      )                               # End capturing ARN for region
      :
      (\\d{12})                       # Account ID. The AWS account ID is a 12-digit number.
      :
      (                               # Begin capturing ARN for "resourcetype/resource"
        (job-queue)                   # Resource type of AWS Batch Job queue
        /                             # Separator between resource type and resource name
        ([\\w-]{1,128})               # Resource name of Job queue can only contain alphanumeric characters, dashes, and underscores. It also must be up to 128 characters long.
      )                               # End capturing ARN for "resourcetype/resource"
    """.trim.r
}

object ArnValidation {
  def apply(key: String): ArnValidation = new ArnValidation(key)
}

class ArnValidation(override val key: String) extends StringRuntimeAttributesValidation(key) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = {
    case WomString(s) => validateArn(s)
  }

  private def validateArn(possibleArn: String): ErrorOr[String] = {
    possibleArn match {
      case arnRegex(_@_*) => possibleArn.validNel
      case _ => "ARN has invalid format".invalidNel
    }
  }

  // Possible ARN formats can be found here
  // https://docs.aws.amazon.com/en_us/general/latest/gr/aws-arns-and-namespaces.html
  // This is quite vague regex, but it allows to support a lot of ARN formats
  protected val arnRegex =
  s"""
      (?x)                            # Turn on comments and whitespace insensitivity
      (arn)                           # Every ARN starts with "arn"
      :
      (                               # Begin capturing ARN for partition
        aws                           # Required part of a partition
        (-[a-z]+){0,2}                # Optional part like "-cn" or "-us-gov". Therefore, the whole partition may look like "aws", "aws-cn", "aws-us-gov"
      )                               # End capturing ARN for partition
      :
      ([a-z0-9-]+)                    # Service name, e.g. "batch", "s3", "aws-marketplace"
      :
      (                               # Begin capturing ARN for region. Not required for some services, e.g. S3
        (
          [a-z]{2}                    # ALPHA-2 county code
          (-[a-z]+){1,2}              # Specific region like "-west" or "-southeast". In case of AWS GovCloud it may be like "-gov-east"
          -\\d                        # Sequence number of the region
        )
        |
        (\\*)                         # Some services like SWF support replacing region with wildcard
      )?                              # End capturing ARN for region.
      :
      (\\d{12})?                      # Account ID. The AWS account ID is a 12-digit number. Not required for some services, e.g. S3
      [:/]{1,2}                       # Separator in this place may not be a simple colon
      (                               # Begin capturing ARN for "resourcetype/resource:(/)qualifier"
        (                             # Begin capturing ARN for resourcetype
          [\\w-]+                     # resourcetype
          [:/]                        # Separator between resourcetype and resource
        )?                            # End capturing ARN for resourcetype
        (                             # Begin capturing ARN for resource
          ([\\w-\\ /.:])+             # Resource name. May have a complex structure like "8ca:auto/my:policy/my"
          (/\\*)?                     # Wildcard that supported by some resources like S2
        )                             # End capturing ARN for resource
        (                             # Begin capturing ARN for optional qualifier
          [:/]                        # Qualifier may be separated both using colon and slash
          [$$\\w]+                    # Qualifier itself
        )?                            # End capturing ARN qualifier. Qualifier is optional
      )                               # End capturing ARN for "resourcetype/resource/qualifier"
    """.trim.r
}

object ZonesValidation extends RuntimeAttributesValidation[Vector[String]] {
  override def key: String = AwsBatchRuntimeAttributes.ZonesKey

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Vector[String]]] = {
    case WomString(s) => s.split("\\s+").toVector.validNel
    case WomArray(womType, value) if womType.memberType == WomStringType =>
      value.map(_.valueString).toVector.validNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be either a whitespace separated String or an Array[String]"
}

object DisksValidation extends RuntimeAttributesValidation[Seq[AwsBatchVolume]] {
  override def key: String = AwsBatchRuntimeAttributes.DisksKey

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[AwsBatchVolume]]] = {
    case WomString(value) => validateLocalDisks(value.split(",\\s*").toSeq)
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      validateLocalDisks(values.map(_.valueString))
  }

  private def validateLocalDisks(disks: Seq[String]): ErrorOr[Seq[AwsBatchVolume]] = {
    val diskNels: Seq[ErrorOr[AwsBatchVolume]] = disks map validateLocalDisk
    val sequenced: ErrorOr[Seq[AwsBatchVolume]] = sequenceNels(diskNels)
    val defaulted: ErrorOr[Seq[AwsBatchVolume]] = addDefault(sequenced)
    defaulted
  }

  private def validateLocalDisk(disk: String): ErrorOr[AwsBatchVolume] = {
    AwsBatchVolume.parse(disk) match {
      case scala.util.Success(attachedDisk) => attachedDisk.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }
  }

  private def sequenceNels(nels: Seq[ErrorOr[AwsBatchVolume]]): ErrorOr[Seq[AwsBatchVolume]] = {
    val emptyDiskNel: ErrorOr[Vector[AwsBatchVolume]] = Vector.empty[AwsBatchVolume].validNel
    val disksNel: ErrorOr[Vector[AwsBatchVolume]] = nels.foldLeft(emptyDiskNel) {
      (acc, v) => (acc, v) mapN { (a, v) => a :+ v }
    }
    disksNel
  }

  private def addDefault(disksNel: ErrorOr[Seq[AwsBatchVolume]]): ErrorOr[Seq[AwsBatchVolume]] = {
    disksNel map {
      case disks if disks.exists(_.name == AwsBatchWorkingDisk.Name) || disks.exists(_.fsType == "efs") => disks
      case disks => disks :+ AwsBatchWorkingDisk.Default
    }
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}
