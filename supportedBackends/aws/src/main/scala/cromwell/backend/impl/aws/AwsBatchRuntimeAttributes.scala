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
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.impl.aws.io.{AwsBatchVolume, AwsBatchWorkingDisk}
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
import wdl4s.parser.MemoryUnit
import wom.types._
import wom.values._
import com.typesafe.config.{ConfigException, ConfigValueFactory}

import scala.util.matching.Regex
import org.slf4j.{Logger, LoggerFactory}
import wom.RuntimeAttributesKeys.GpuKey

import scala.util.{Failure, Success, Try}
import scala.jdk.CollectionConverters._

/**
 * Attributes that are provided to the job at runtime
 *
 * @param cpu number of vCPU
 * @param gpuCount number of gpu
 * @param zones the aws availability zones to run in
 * @param memory memory to allocate
 * @param disks a sequence of disk volumes
 * @param dockerImage the name of the docker image that the job will run in
 * @param queueArn the arn of the AWS Batch queue that the job will be submitted to
 * @param failOnStderr should the job fail if something is logged to `stderr`
 * @param continueOnReturnCode decides if a job continues on receiving a specific return code
 * @param noAddress is there no address
 * @param scriptS3BucketName the s3 bucket where the execution command or script will be written and, from there, fetched into the container and executed
 * @param fileSystem the filesystem type, default is "s3"
 * @param awsBatchRetryAttempts number of attempts that AWS Batch will retry the task if it fails
 * @param awsBatchEvaluateOnExit Evaluate on exit strategy setting for AWS batch retry
 * @param ulimits ulimit values to be passed to the container
 * @param efsDelocalize should we delocalize efs files to s3
 * @param efsMakeMD5 should we make a sibling md5 file as part of the job
 * @param tagResources should we tag resources
 * @param logGroupName the CloudWatch log group name to write logs to
 * @param additionalTags a map of tags to add to the AWS Batch job submission
  * @param propagateTags should we propagate tags to ECS tasks underlying batch jobs
 */
case class AwsBatchRuntimeAttributes(cpu: Int Refined Positive,
                                     gpuCount: Int,
                                     zones: Vector[String],
                                     memory: MemorySize,
                                     disks: Seq[AwsBatchVolume],
                                     dockerImage: String,
                                     queueArn: String,
                                     failOnStderr: Boolean,
                                     continueOnReturnCode: ContinueOnReturnCode,
                                     noAddress: Boolean,
                                     scriptS3BucketName: String,
                                     awsBatchRetryAttempts: Int,
                                     awsBatchEvaluateOnExit: Vector[Map[String, String]],
                                     ulimits: Vector[Map[String, String]],
                                     efsDelocalize: Boolean,
                                     efsMakeMD5: Boolean,
                                     sharedMemorySize: MemorySize,
                                     jobTimeout: Int,
                                     logGroupName: String,
                                     additionalTags: Map[String, String],
                                     fuseMount: Boolean,
                                     fileSystem: String = "s3",
                                     tagResources: Boolean = false,
                                     propagateTags: Boolean = false
)

object AwsBatchRuntimeAttributes {
  val Log: Logger = LoggerFactory.getLogger(this.getClass)
  val QueueArnKey = "queueArn"

  val scriptS3BucketKey = "scriptBucketName"

  val awsBatchRetryAttemptsKey = "awsBatchRetryAttempts"

  val awsBatchEvaluateOnExitKey = "awsBatchEvaluateOnExit"

  val defaultSharedMemorySize = MemorySize(64, MemoryUnit.MB)

  private val awsBatchEvaluateOnExitDefault =
    WomArray(WomArrayType(WomMapType(WomStringType, WomStringType)), Vector(WomMap(Map.empty[WomValue, WomValue])))

  val awsBatchefsDelocalizeKey = "efsDelocalize"
  val awsBatchefsMakeMD5Key = "efsMakeMD5"
  val tagResourcesKey = "tagResources"
  val propagateResourcesKey = "propagateTags"
  val ZonesKey = "zones"
  private val ZonesDefaultValue = WomString("us-east-1a")

  val NoAddressKey = "noAddress"
  private val noAddressValidationInstance = new BooleanRuntimeAttributesValidation(NoAddressKey)
  private val NoAddressDefaultValue = WomBoolean(false)

  // TODO: Determine good volume format
  val DisksKey = "disks"
  private val DisksDefaultValue = WomString(s"${AwsBatchWorkingDisk.Name}")

  private val MemoryDefaultValue = "2 GB"

  private val logGroupNameKey = "logGroupName"
  private val logGroupNameValidationInstance = new StringRuntimeAttributesValidation(logGroupNameKey)
  private val LogGroupNameDefaultValue = WomString("/aws/batch/job")

  private val additionalTagsKey = "additionalTags"

  val UlimitsKey = "ulimits"
  private val jobTimeoutKey = "jobTimeout"

  val fuseMountKey = "fuseMount"
  private val fuseMountValidationInstance = new BooleanRuntimeAttributesValidation(fuseMountKey)
  private val fuseMountDefaultValue = WomBoolean(false)

  private val UlimitsDefaultValue =
    WomArray(WomArrayType(WomMapType(WomStringType, WomStringType)), Vector(WomMap(Map.empty[WomValue, WomValue])))

  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] =
    CpuValidation.instance
      .withDefault(CpuValidation.configDefaultWomValue(runtimeConfig) getOrElse CpuValidation.defaultMin)

  private def gpuCountValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] =
    PosIntValidation(GpuKey).withDefault(
      PosIntValidation(GpuKey).configDefaultWomValue(runtimeConfig).getOrElse(WomInteger(0))
    )

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) =
    ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def disksValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Seq[AwsBatchVolume]] =
    DisksValidation
      .withDefault(DisksValidation.configDefaultWomValue(runtimeConfig) getOrElse DisksDefaultValue)

  private def zonesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Vector[String]] =
    ZonesValidation
      .withDefault(ZonesValidation.configDefaultWomValue(runtimeConfig) getOrElse ZonesDefaultValue)

  private def memoryValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[MemorySize] =
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.MemoryKey,
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, runtimeConfig) getOrElse MemoryDefaultValue
    )

  private def noAddressValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    noAddressValidationInstance
      .withDefault(noAddressValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse NoAddressDefaultValue)

  private def logGroupNameValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[String] =
    logGroupNameValidationInstance
      .withDefault(
        logGroupNameValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse LogGroupNameDefaultValue
      )

  private def sharedMemorySizeValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[MemorySize] =
    MemoryValidation.withDefaultMemory(
      RuntimeAttributesKeys.sharedMemoryKey,
      MemoryValidation.configDefaultString(RuntimeAttributesKeys.sharedMemoryKey,
                                           runtimeConfig
      ) getOrElse defaultSharedMemorySize.toString
    )

  private def fuseMountValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    fuseMountValidationInstance.withDefault(
      fuseMountValidationInstance.configDefaultWomValue(runtimeConfig) getOrElse fuseMountDefaultValue
    )

  private def jobTimeoutValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] =
    PosIntValidation(jobTimeoutKey, minValue = 60).withDefault(
      PosIntValidation(jobTimeoutKey, minValue = 60).configDefaultWomValue(runtimeConfig).getOrElse(WomInteger(0))
    )

  private def scriptS3BucketNameValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[String] =
    ScriptS3BucketNameValidation(scriptS3BucketKey).withDefault(
      ScriptS3BucketNameValidation(scriptS3BucketKey)
        .configDefaultWomValue(runtimeConfig)
        .getOrElse(throw new RuntimeException("scriptBucketName is required"))
    )

  // As of WDL 1.1 these two are aliases of each other
  private val dockerValidation: OptionalRuntimeAttributesValidation[Containers] = DockerValidation.instance
  private val containerValidation: OptionalRuntimeAttributesValidation[Containers] = ContainerValidation.instance

  private def queueArnValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[String] =
    QueueArnValidation.withDefault(
      QueueArnValidation.configDefaultWomValue(runtimeConfig) getOrElse
        (throw new RuntimeException("queueArn is required"))
    )

  private def awsBatchRetryAttemptsValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int] =
    AwsBatchRetryAttemptsValidation(awsBatchRetryAttemptsKey).withDefault(
      AwsBatchRetryAttemptsValidation(awsBatchRetryAttemptsKey)
        .configDefaultWomValue(runtimeConfig)
        .getOrElse(WomInteger(0))
    )

  def awsBatchEvaluateOnExitValidation(
    runtimeConfig: Option[Config]
  ): RuntimeAttributesValidation[Vector[Map[String, String]]] =
    AwsBatchEvaluateOnExitValidation
      .withDefault(AwsBatchEvaluateOnExitValidation.fromConfig(runtimeConfig).getOrElse(awsBatchEvaluateOnExitDefault))

  private def awsBatchefsDelocalizeValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    AwsBatchefsDelocalizeValidation(awsBatchefsDelocalizeKey).withDefault(
      AwsBatchefsDelocalizeValidation(awsBatchefsDelocalizeKey)
        .configDefaultWomValue(runtimeConfig)
        .getOrElse(WomBoolean(false))
    )

  private def awsBatchefsMakeMD5Validation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    AwsBatchefsMakeMD5Validation(awsBatchefsMakeMD5Key).withDefault(
      AwsBatchefsMakeMD5Validation(awsBatchefsMakeMD5Key)
        .configDefaultWomValue(runtimeConfig)
        .getOrElse(WomBoolean(false))
    )

  private def awsBatchtagResourcesValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    AwsBatchtagResourcesValidation(AwsBatchRuntimeAttributes.tagResourcesKey).withDefault(
      AwsBatchtagResourcesValidation(AwsBatchRuntimeAttributes.tagResourcesKey)
        .configDefaultWomValue(runtimeConfig)
        .getOrElse(WomBoolean(false))
    )

  private def awsBatchPropagateTagsValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] =
    AwsBatchtagResourcesValidation(AwsBatchRuntimeAttributes.propagateResourcesKey).withDefault(
      AwsBatchtagResourcesValidation(AwsBatchRuntimeAttributes.propagateResourcesKey)
        .configDefaultWomValue(runtimeConfig)
        .getOrElse(WomBoolean(false))
    )

  private def ulimitsValidation(
    runtimeConfig: Option[Config]
  ): RuntimeAttributesValidation[Vector[Map[String, String]]] =
    UlimitsValidation.withDefault(UlimitsValidation.configDefaultWomValue(runtimeConfig) getOrElse UlimitsDefaultValue)

  // routine that aggregates disks from default-runtime-attributes and efs.
  def aggregateDisksInRuntimeConfig(configuration: AwsBatchConfiguration): Option[Config] = {
    // get the runtime attributes out of full config
    val rtc: Config = configuration.runtimeConfig match {
      case Some(x) => x
      case None =>
        Log.error("aws-batch default runtime attributes not found.")
        throw new RuntimeException(
          "Default runtime attributes not found in config. This should always be present for AWS ! "
        )
    }
    // efs disk configs is set (can be None) in configuration
    val efs_disks = configuration.efsMntPoint.getOrElse("")
    // additional disks optionally set as runtime attributes
    val rtc_disks =
      try
        rtc
          .getAnyRef(AwsBatchRuntimeAttributes.DisksKey)
          .asInstanceOf[String] // just to prevent complaints about var/val
      catch {
        case _: ConfigException.Missing =>
          "local-disk"
      }
    // combine and remove empty values and remove empty values
    val disks = s"${efs_disks},${rtc_disks}".split(",").toSet.filterNot(_.isEmpty).filterNot(_.isEmpty).mkString(",")
    val runtimeConfig = Some(rtc.withValue(AwsBatchRuntimeAttributes.DisksKey, ConfigValueFactory.fromAnyRef(disks)))
    return runtimeConfig
  }

  def runtimeAttributesBuilder(configuration: AwsBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = aggregateDisksInRuntimeConfig(configuration)
    def validationsS3backend = StandardValidatedRuntimeAttributesBuilder
      .default(runtimeConfig)
      .withValidation(
        cpuValidation(runtimeConfig),
        gpuCountValidation(runtimeConfig),
        disksValidation(runtimeConfig),
        zonesValidation(runtimeConfig),
        memoryValidation(runtimeConfig),
        noAddressValidation(runtimeConfig),
        dockerValidation,
        containerValidation,
        queueArnValidation(runtimeConfig),
        scriptS3BucketNameValidation(runtimeConfig),
        logGroupNameValidation(runtimeConfig),
        awsBatchRetryAttemptsValidation(runtimeConfig),
        awsBatchEvaluateOnExitValidation(runtimeConfig),
        ulimitsValidation(runtimeConfig),
        awsBatchefsDelocalizeValidation(runtimeConfig),
        awsBatchefsMakeMD5Validation(runtimeConfig),
        awsBatchtagResourcesValidation(runtimeConfig),
        awsBatchPropagateTagsValidation(runtimeConfig),
        sharedMemorySizeValidation(runtimeConfig),
        fuseMountValidation(runtimeConfig),
        jobTimeoutValidation(runtimeConfig)
      )
    def validationsLocalBackend = StandardValidatedRuntimeAttributesBuilder
      .default(runtimeConfig)
      .withValidation(
        cpuValidation(runtimeConfig),
        gpuCountValidation(runtimeConfig),
        disksValidation(runtimeConfig),
        zonesValidation(runtimeConfig),
        memoryValidation(runtimeConfig),
        noAddressValidation(runtimeConfig),
        dockerValidation,
        containerValidation,
        queueArnValidation(runtimeConfig),
        logGroupNameValidation(runtimeConfig),
        awsBatchRetryAttemptsValidation(runtimeConfig),
        awsBatchEvaluateOnExitValidation(runtimeConfig),
        ulimitsValidation(runtimeConfig),
        awsBatchefsDelocalizeValidation(runtimeConfig),
        awsBatchefsMakeMD5Validation(runtimeConfig),
        awsBatchtagResourcesValidation(runtimeConfig),
        awsBatchPropagateTagsValidation(runtimeConfig),
        sharedMemorySizeValidation(runtimeConfig),
        fuseMountValidation(runtimeConfig),
        jobTimeoutValidation(runtimeConfig)
      )

    configuration.fileSystem match {
      case AWSBatchStorageSystems.s3 => validationsS3backend

      case _ => validationsLocalBackend
    }
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes,
            runtimeAttrsConfig: Option[Config],
            fileSystem: String
  ): AwsBatchRuntimeAttributes = {
    val cpu: Int Refined Positive =
      RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val gpuCount: Int =
      RuntimeAttributesValidation.extract(gpuCountValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val zones: Vector[String] = RuntimeAttributesValidation.extract(ZonesValidation, validatedRuntimeAttributes)
    val memory: MemorySize =
      RuntimeAttributesValidation.extract(memoryValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val disks: Seq[AwsBatchVolume] =
      RuntimeAttributesValidation.extract(disksValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val docker: String = Containers.extractContainer(validatedRuntimeAttributes)
    val queueArn: String =
      RuntimeAttributesValidation.extract(queueArnValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(
      continueOnReturnCodeValidation(runtimeAttrsConfig),
      validatedRuntimeAttributes
    )
    val noAddress: Boolean =
      RuntimeAttributesValidation.extract(noAddressValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val scriptS3BucketName = fileSystem match {
      case AWSBatchStorageSystems.s3 =>
        RuntimeAttributesValidation.extract(scriptS3BucketNameValidation(runtimeAttrsConfig),
                                            validatedRuntimeAttributes
        )
      case _ => ""
    }
    val logGroupName: String =
      RuntimeAttributesValidation.extract(logGroupNameValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val additionalTags: Map[String, String] = runtimeAttrsConfig
      .collect {
        case config if config.hasPath(additionalTagsKey) =>
          config
            .getObject(additionalTagsKey)
            .entrySet()
            .asScala
            .map(e => e.getKey -> e.getValue.unwrapped().toString)
            .toMap
      }
      .getOrElse(Map.empty[String, String])

    val awsBatchRetryAttempts: Int = RuntimeAttributesValidation.extract(
      awsBatchRetryAttemptsValidation(runtimeAttrsConfig),
      validatedRuntimeAttributes
    )
    val awsBatchEvaluateOnExit: Vector[Map[String, String]] = RuntimeAttributesValidation.extract(
      awsBatchEvaluateOnExitValidation(runtimeAttrsConfig),
      validatedRuntimeAttributes
    )

    val ulimits: Vector[Map[String, String]] =
      RuntimeAttributesValidation.extract(ulimitsValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val efsDelocalize: Boolean = RuntimeAttributesValidation.extract(
      awsBatchefsDelocalizeValidation(runtimeAttrsConfig),
      validatedRuntimeAttributes
    )
    val efsMakeMD5: Boolean =
      RuntimeAttributesValidation.extract(awsBatchefsMakeMD5Validation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val tagResources: Boolean = RuntimeAttributesValidation.extract(awsBatchtagResourcesValidation(runtimeAttrsConfig),
                                                                    validatedRuntimeAttributes
    )
    val propagateTags: Boolean = RuntimeAttributesValidation.extract(
      awsBatchPropagateTagsValidation(runtimeAttrsConfig),
      validatedRuntimeAttributes
    )
    val sharedMemorySize: MemorySize =
      RuntimeAttributesValidation.extract(sharedMemorySizeValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val jobTimeout: Int =
      RuntimeAttributesValidation.extract(jobTimeoutValidation(runtimeAttrsConfig), validatedRuntimeAttributes)
    val fuseMount: Boolean =
      RuntimeAttributesValidation.extract(fuseMountValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    new AwsBatchRuntimeAttributes(
      cpu,
      gpuCount,
      zones,
      memory,
      disks,
      docker,
      queueArn,
      failOnStderr,
      continueOnReturnCode,
      noAddress,
      scriptS3BucketName,
      awsBatchRetryAttempts,
      awsBatchEvaluateOnExit,
      ulimits,
      efsDelocalize,
      efsMakeMD5,
      sharedMemorySize,
      jobTimeout,
      logGroupName,
      additionalTags,
      fuseMount,
      fileSystem,
      tagResources,
      propagateTags
    )
  }
}

object ScriptS3BucketNameValidation {
  def apply(key: String): ScriptS3BucketNameValidation = new ScriptS3BucketNameValidation(key)
}

class ScriptS3BucketNameValidation(key: String) extends StringRuntimeAttributesValidation(key) {

  // a reasonable but not perfect regex for a bucket. see https://stackoverflow.com/a/50484916/3573553
  protected val s3BucketNameRegex: Regex =
    "(?=^.{3,63}$)(?!^(\\d+\\.)+\\d+$)(^(([a-z0-9]|[a-z0-9][a-z0-9\\-]*[a-z0-9])\\.)*([a-z0-9]|[a-z0-9][a-z0-9\\-]*[a-z0-9])$)".r

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = { case WomString(s) =>
    validateBucketName(s)
  }

  private def validateBucketName(possibleBucketName: String): ErrorOr[String] =
    possibleBucketName match {
      case s3BucketNameRegex(_ @_*) => possibleBucketName.validNel
      case _ => "The Script Bucket name has an invalid s3 bucket format".invalidNel
    }
}

object QueueArnValidation extends ArnValidation(AwsBatchRuntimeAttributes.QueueArnKey) {
  // queue arn format can be found here
  // https://docs.aws.amazon.com/en_us/general/latest/gr/aws-arns-and-namespaces.html#arn-syntax-batch
  // arn:aws:batch:region:account-id:job-queue/queue-name
  override protected val arnRegex: Regex =
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
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[String]] = { case WomString(s) =>
    validateArn(s)
  }

  private def validateArn(possibleArn: String): ErrorOr[String] =
    possibleArn match {
      case arnRegex(_ @_*) => possibleArn.validNel
      case _ => "ARN has invalid format".invalidNel
    }

  // Possible ARN formats can be found here
  // https://docs.aws.amazon.com/en_us/general/latest/gr/aws-arns-and-namespaces.html
  // This is quite vague regex, but it allows to support a lot of ARN formats
  protected val arnRegex: Regex =
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

  override def coercion: Iterable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

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

  override def coercion: Iterable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

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

  private def validateLocalDisk(disk: String): ErrorOr[AwsBatchVolume] =
    AwsBatchVolume.parse(disk) match {
      case scala.util.Success(attachedDisk) => attachedDisk.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }

  private def sequenceNels(nels: Seq[ErrorOr[AwsBatchVolume]]): ErrorOr[Seq[AwsBatchVolume]] = {
    val emptyDiskNel: ErrorOr[Vector[AwsBatchVolume]] = Vector.empty[AwsBatchVolume].validNel
    val disksNel: ErrorOr[Vector[AwsBatchVolume]] = nels.foldLeft(emptyDiskNel) { (acc, v) =>
      (acc, v) mapN { (a, v) => a :+ v }
    }
    disksNel
  }

  private def addDefault(disksNel: ErrorOr[Seq[AwsBatchVolume]]): ErrorOr[Seq[AwsBatchVolume]] =
    disksNel map {
      case disks if disks.exists(_.name == AwsBatchWorkingDisk.Name) || disks.exists(_.fsType == "efs") => disks
      case disks => disks :+ AwsBatchWorkingDisk.Default
    }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}

object PosIntValidation {
  def apply(key: String, minValue: Int = 0): PosIntValidation = new PosIntValidation(key, minValue)
}

class PosIntValidation(key: String, minValue: Int = 0) extends IntRuntimeAttributesValidation(key) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Int]] = {
    case womValue if WomIntegerType.coerceRawValue(womValue).isSuccess =>
      WomIntegerType.coerceRawValue(womValue).get match {
        case WomInteger(value) =>
          // allow the default value of zero.
          if (value.toInt < minValue && value.toInt != 0)
            s"Expecting $key runtime attribute value greater than or equal to $minValue".invalidNel
          else
            value.toInt.validNel
      }
  }
}

object AwsBatchRetryAttemptsValidation {
  def apply(key: String): AwsBatchRetryAttemptsValidation = new AwsBatchRetryAttemptsValidation(key)
}

class AwsBatchRetryAttemptsValidation(key: String) extends IntRuntimeAttributesValidation(key) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Int]] = {
    case womValue if WomIntegerType.coerceRawValue(womValue).isSuccess =>
      WomIntegerType.coerceRawValue(womValue).get match {
        case WomInteger(value) =>
          if (value.toInt < 0)
            s"Expecting $key runtime attribute value greater than or equal to 0".invalidNel
          else if (value.toInt > 10)
            s"Expecting $key runtime attribute value lower than or equal to 10".invalidNel
          else
            value.toInt.validNel
      }
  }

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be an Integer"
}

object AwsBatchEvaluateOnExitValidation extends RuntimeAttributesValidation[Vector[Map[String, String]]] {

  val requiredKey = "action"
  private val acceptedKeys = Set(requiredKey, "onExitCode", "onReason", "onStatusReason")

  def fromConfig(runtimeConfig: Option[Config]): Option[WomValue] = {
    val config = runtimeConfig match {
      case Some(value) =>
        Try(value.getObjectList(key)) match {
          case Failure(_) => None
          case Success(value) =>
            Some(value.asScala.map {
              _.unwrapped().asScala.toMap
            }.toList)
        }
      case _ => None
    }

    config match {
      case Some(value) =>
        Some(AwsBatchEvaluateOnExitValidation.coercion collectFirst {
          case womType if womType.coerceRawValue(value).isSuccess => womType.coerceRawValue(value).get
        } getOrElse {
          BadDefaultAttribute(WomString(value.toString))
        })
      case None => None
    }
  }

  override def coercion: Iterable[WomType] =
    Set(WomStringType, WomArrayType(WomMapType(WomStringType, WomStringType)))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Vector[Map[String, String]]]] = {
    case WomArray(womType, value) if womType.memberType == WomMapType(WomStringType, WomStringType) =>
      check_maps(value.toVector)
    case WomMap(_, _) => "!!! ERROR1".invalidNel
  }

  private def check_maps(
    maps: Vector[WomValue]
  ): ErrorOr[Vector[Map[String, String]]] = {
    val entryNels: Vector[ErrorOr[Map[String, String]]] = maps.map {
      case WomMap(_, value) => check_keys(value)
      case _ => "!!! ERROR2".invalidNel
    }
    val sequenced: ErrorOr[Vector[Map[String, String]]] = sequenceNels(
      entryNels
    )
    sequenced
  }

  private def validateActionKey(dict: Map[WomValue, WomValue]): ErrorOr[Map[String, String]] = {
    val validCondition = Set("retry", "exit")
    val convertedMap = dict
      .map { case (WomString(k), WomString(v)) =>
        (k, v)
      // case _ => "!!! ERROR3".invalidNel
      }
    if (
      convertedMap.exists { case (key, value) =>
        key.toLowerCase == requiredKey && validCondition.contains(value.toLowerCase)
      }
    ) {
      convertedMap.validNel
    } else {
      s"Missing or invalid $requiredKey key/value for runtime attribute: $key. Refer to https://docs.aws.amazon.com/batch/latest/APIReference/API_RetryStrategy.html".invalidNel
    }
  }

  private def check_keys(
    dict: Map[WomValue, WomValue]
  ): ErrorOr[Map[String, String]] = {
    val map_keys = dict.keySet.map(_.valueString.toLowerCase)
    val unrecognizedKeys =
      map_keys.diff(acceptedKeys.map(x => x.toLowerCase))
    if (!dict.nonEmpty) {
      Map.empty[String, String].validNel
    } else if (unrecognizedKeys.nonEmpty) {
      s"Invalid keys in $key runtime attribute: $unrecognizedKeys. Only $acceptedKeys are accepted. Refer to https://docs.aws.amazon.com/batch/latest/APIReference/API_RetryStrategy.html".invalidNel
    } else {
      validateActionKey(dict)
    }
  }

  private def sequenceNels(
    nels: Vector[ErrorOr[Map[String, String]]]
  ): ErrorOr[Vector[Map[String, String]]] = {
    val emptyNel: ErrorOr[Vector[Map[String, String]]] =
      Vector.empty[Map[String, String]].validNel
    val seqNel: ErrorOr[Vector[Map[String, String]]] =
      nels.foldLeft(emptyNel) { (acc, v) =>
        (acc, v) mapN { (a, v) => a :+ v }
      }
    seqNel
  }

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be defined"

  /**
   * Returns the key of the runtime attribute.
   *
   * @return The key of the runtime attribute.
   */
  override def key: String = AwsBatchRuntimeAttributes.awsBatchEvaluateOnExitKey
}

object AwsBatchefsDelocalizeValidation {
  def apply(key: String): AwsBatchefsDelocalizeValidation = new AwsBatchefsDelocalizeValidation(key)
}

class AwsBatchefsDelocalizeValidation(key: String) extends BooleanRuntimeAttributesValidation(key) {

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be an Boolean"
}

object AwsBatchefsMakeMD5Validation {
  def apply(key: String): AwsBatchefsMakeMD5Validation = new AwsBatchefsMakeMD5Validation(key)
}

class AwsBatchefsMakeMD5Validation(key: String) extends BooleanRuntimeAttributesValidation(key) {

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be an Boolean"
}

object AwsBatchtagResourcesValidation {
  def apply(key: String): AwsBatchtagResourcesValidation = new AwsBatchtagResourcesValidation(key)
}

class AwsBatchtagResourcesValidation(key: String) extends BooleanRuntimeAttributesValidation(key) {

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be a Boolean"
}

object UlimitsValidation extends RuntimeAttributesValidation[Vector[Map[String, String]]] {
  override def key: String = AwsBatchRuntimeAttributes.UlimitsKey

  override def coercion: Iterable[WomType] =
    Set(WomStringType, WomArrayType(WomMapType(WomStringType, WomStringType)))

  var accepted_keys = Set("name", "softLimit", "hardLimit")

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Vector[Map[String, String]]]] = {
    case WomArray(womType, value) if womType.memberType == WomMapType(WomStringType, WomStringType) =>
      check_maps(value.toVector)
    case WomMap(_, _) => "!!! ERROR1".invalidNel

  }

  private def check_maps(
    maps: Vector[WomValue]
  ): ErrorOr[Vector[Map[String, String]]] = {
    val entryNels: Vector[ErrorOr[Map[String, String]]] = maps.map {
      case WomMap(_, value) => check_keys(value)
      case _ => "!!! ERROR2".invalidNel
    }
    val sequenced: ErrorOr[Vector[Map[String, String]]] = sequenceNels(
      entryNels
    )
    sequenced
  }

  private def check_keys(
    dict: Map[WomValue, WomValue]
  ): ErrorOr[Map[String, String]] = {
    val map_keys = dict.keySet.map(_.valueString).toSet
    val unrecognizedKeys =
      accepted_keys.diff(map_keys) union map_keys.diff(accepted_keys)

    if (!dict.nonEmpty) {
      Map.empty[String, String].validNel
    } else if (unrecognizedKeys.nonEmpty) {
      s"Invalid keys in $key runtime attribute. Refer to 'ulimits' section on https://docs.aws.amazon.com/batch/latest/userguide/job_definition_parameters.html#containerProperties".invalidNel
    } else {
      dict
        .collect { case (WomString(k), WomString(v)) =>
          (k, v)
        // case _ => "!!! ERROR3".invalidNel
        }
        .toMap
        .validNel
    }
  }

  private def sequenceNels(
    nels: Vector[ErrorOr[Map[String, String]]]
  ): ErrorOr[Vector[Map[String, String]]] = {
    val emptyNel: ErrorOr[Vector[Map[String, String]]] =
      Vector.empty[Map[String, String]].validNel
    val seqNel: ErrorOr[Vector[Map[String, String]]] =
      nels.foldLeft(emptyNel) { (acc, v) =>
        (acc, v) mapN { (a, v) => a :+ v }
      }
    seqNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be an Array[Map[String, String]]"
}
