package cromwell.backend.impl.bcs

import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import com.typesafe.config.Config
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.types.{WdlArrayType, WdlStringType, WdlType}
import wdl4s.wdl.values.{WdlArray, WdlBoolean, WdlInteger, WdlString, WdlValue}

import scala.util.{Failure, Success, Try}

trait OptionalWithDefault[A] {
  this: RuntimeAttributesValidation[A] =>
  protected val config: Option[Config]

  override protected def staticDefaultOption: Option[WdlValue] = {
    Try(this.configDefaultWdlValue(config)) match {
      case Success(value: Option[WdlValue]) => value
      case Failure(_) => None
    }
  }
}

case class BcsRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                dockerImage: Option[String],
                                dockerPath: Option[String],
                                failOnStderr: Boolean,
                                mounts: Option[Seq[BcsMount]],
                                userData: Option[Seq[BcsUserData]],
                                clusterId: Option[String],
                                resourceType: Option[String],
                                instanceType: Option[String],
                                imageId: Option[String],
                                reserveOnFail: Option[Boolean],
                                autoReleaseJob: Option[Boolean],
                                workerPath: Option[String],
                                timeout: Option[Int],
                                verbose: Option[Boolean])

object BcsRuntimeAttributes {

  val MountsKey = "mounts"
  val UserDataKey = "userData"
  val MountsDefaultValue = WdlString("")
  val ReserveOnFailKey = "reserveOnFail"
  val ReserveOnFailDefault = false
  val AutoReleaseJobKey = "autoReleaseJob"
  val AutoReleaseJobDefault = WdlBoolean(true)
  val TimeoutKey = "timeout"
  val TimeoutDefault = WdlInteger(21600)
  val VerboseKey = "verbose"

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def dockerImageValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = DockerImageValidation.optionalWithDefault(runtimeConfig)

  private def dockerPathValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = DockerPathValidation.optionalWithDefault(runtimeConfig)


  private def clusterIdValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = ClusterIdValidation.optionalWithDefault(runtimeConfig)

  private def resourceTypeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = ResourceTypeValidation.optionalWithDefault(runtimeConfig)

  private def instanceTypeValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = InstanceTypeValidation.optionalWithDefault(runtimeConfig)

  private def imageIdValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = ImageIdValidation.optionalWithDefault(runtimeConfig)

  private def userDataValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsUserData]] = UserDataValidation.optionalWithDefault(runtimeConfig)

  private def reserveOnFailValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = ReserveOnFailValidation.optionalWithDefault(runtimeConfig)


  private def autoReleaseJobValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = AutoReleaseJobValidation.optionalWithDefault(runtimeConfig)

  private def mountsValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsMount]] = MountsValidation.optionalWithDefault(runtimeConfig)

  private def workerPathValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = WorkerPathValidation.optionalWithDefault(runtimeConfig)

  private def timeoutValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int] = TimeoutValidation.optionalWithDefault(runtimeConfig)

  private def verboseValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = VerboseValidation.optionalWithDefault(runtimeConfig)

  def runtimeAttributesBuilder(backendRuntimeConfig: Option[Config]): StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default(backendRuntimeConfig).withValidation(
      mountsValidation(backendRuntimeConfig),
      userDataValidation(backendRuntimeConfig),
      clusterIdValidation(backendRuntimeConfig),
      dockerImageValidation(backendRuntimeConfig),
      dockerPathValidation(backendRuntimeConfig),
      resourceTypeValidation(backendRuntimeConfig),
      instanceTypeValidation(backendRuntimeConfig),
      imageIdValidation(backendRuntimeConfig),
      reserveOnFailValidation(backendRuntimeConfig),
      autoReleaseJobValidation(backendRuntimeConfig),
      workerPathValidation(backendRuntimeConfig),
      timeoutValidation(backendRuntimeConfig),
      verboseValidation(backendRuntimeConfig)
    )

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, backendRuntimeConfig: Option[Config]): BcsRuntimeAttributes = {
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val mounts: Option[Seq[BcsMount]] = RuntimeAttributesValidation.extractOption(mountsValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val userData: Option[Seq[BcsUserData]] = RuntimeAttributesValidation.extractOption(userDataValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val clusterId: Option[String] = RuntimeAttributesValidation.extractOption(clusterIdValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val dockerImage: Option[String] = RuntimeAttributesValidation.extractOption(dockerImageValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val dockerPath: Option[String] = RuntimeAttributesValidation.extractOption(dockerPathValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val resourceType: Option[String] = RuntimeAttributesValidation.extractOption(resourceTypeValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val instanceType: Option[String] = RuntimeAttributesValidation.extractOption(instanceTypeValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val imageId: Option[String] = RuntimeAttributesValidation.extractOption(imageIdValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val reserveOnFail: Option[Boolean] = RuntimeAttributesValidation.extractOption(reserveOnFailValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val autoReleaseJob: Option[Boolean] = RuntimeAttributesValidation.extractOption(autoReleaseJobValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val workerPath: Option[String] = RuntimeAttributesValidation.extractOption(workerPathValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val timeout: Option[Int] = RuntimeAttributesValidation.extractOption(timeoutValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val verbose: Option[Boolean] = RuntimeAttributesValidation.extractOption(verboseValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    new BcsRuntimeAttributes(
      continueOnReturnCode,
      dockerImage,
      dockerPath,
      failOnStderr,
      mounts,
      userData,
      clusterId,
      resourceType,
      instanceType,
      imageId,
      reserveOnFail,
      autoReleaseJob,
      workerPath,
      timeout,
      verbose
    )
  }
}

object MountsValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsMount]] = new MountsValidation(config).optional
}

class MountsValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[Seq[BcsMount]] with OptionalWithDefault[Seq[BcsMount]] {
  override def key: String = BcsRuntimeAttributes.MountsKey

  override def coercion: Traversable[WdlType] = Set(WdlStringType, WdlArrayType(WdlStringType))

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[Seq[BcsMount]]] = {
    case WdlString(value) => validateMounts(value.split(",\\s*").toSeq)
    case WdlArray(wdlType, values) if wdlType.memberType == WdlStringType =>
      validateMounts(values.map(_.valueString))
  }

  private def validateMounts(mounts: Seq[String]): ErrorOr[Seq[BcsMount]] = {
    val mountNels: Seq[ErrorOr[BcsMount]] = mounts filter { s => !s.trim().isEmpty } map validateMounts
    val sequenced: ErrorOr[Seq[BcsMount]] = sequenceNels(mountNels)
    sequenced
  }

  private def validateMounts(mount: String): ErrorOr[BcsMount] = {
    BcsMount.parse(mount) match {
      case scala.util.Success(mnt) => mnt.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }
  }

  private def sequenceNels(nels: Seq[ErrorOr[BcsMount]]): ErrorOr[Seq[BcsMount]] = {
    val emptyMountNel = Vector.empty[BcsMount].validNel[String]
    val mountsNel: ErrorOr[Vector[BcsMount]] = nels.foldLeft(emptyMountNel) {
      (acc, v) => (acc |@| v) map { (a, v) => a :+ v }
    }
    mountsNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}

object UserDataValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsUserData]] = new UserDataValidation(config).optional
}

class UserDataValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[Seq[BcsUserData]] with OptionalWithDefault[Seq[BcsUserData]] {
  override def key: String = BcsRuntimeAttributes.UserDataKey

  override def coercion: Traversable[WdlType] = Set(WdlStringType, WdlArrayType(WdlStringType))

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[Seq[BcsUserData]]] = {
    case WdlString(value) => validateUserData(value.split(",\\s*").toSeq)
    case WdlArray(wdlType, values) if wdlType.memberType == WdlStringType =>
      validateUserData(values.map(_.valueString))
  }

  private def validateUserData(mounts: Seq[String]): ErrorOr[Seq[BcsUserData]] = {
    val userDataNels: Seq[ErrorOr[BcsUserData]] = mounts filter { s => !s.trim().isEmpty } map validateUserData
    val sequenced: ErrorOr[Seq[BcsUserData]] = sequenceNels(userDataNels)
    sequenced
  }

  private def validateUserData(data: String): ErrorOr[BcsUserData] = {
    BcsUserData.parse(data) match {
      case scala.util.Success(userData) => userData.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }
  }

  private def sequenceNels(nels: Seq[ErrorOr[BcsUserData]]): ErrorOr[Seq[BcsUserData]] = {
    val emptyDataNel = Vector.empty[BcsUserData].validNel[String]
    val datasNel: ErrorOr[Vector[BcsUserData]] = nels.foldLeft(emptyDataNel) {
      (acc, v) => (acc |@| v) map { (a, v) => a :+ v }
    }
    datasNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}

object ClusterIdValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new ClusterIdValidation(config).optional
}

class ClusterIdValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("cluster") with OptionalWithDefault[String] {
  override protected def usedInCallCaching: Boolean = false

  override protected def missingValueMessage: String = "Can't find an attribute value for key cluster-id"

  override protected def invalidValueMessage(value: WdlValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object DockerImageValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new DockerImageValidation(config).optional
}

class DockerImageValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("dockerImage") with OptionalWithDefault[String] {
  override protected def usedInCallCaching: Boolean = false

  override protected def missingValueMessage: String = "Can't find an attribute value for key docker image"

  override protected def invalidValueMessage(value: WdlValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object DockerPathValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new DockerPathValidation(config).optional
}

class DockerPathValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("dockerPath") with OptionalWithDefault[String] {
  override protected def usedInCallCaching: Boolean = false

  override protected def missingValueMessage: String = "Can't find an attribute value for key docker path"

  override protected def invalidValueMessage(value: WdlValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object ResourceTypeValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new ResourceTypeValidation(config).optional
}

class ResourceTypeValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("resourceType") with OptionalWithDefault[String] {
  override protected def usedInCallCaching: Boolean = false

  override protected def missingValueMessage: String = "Can't find an attribute value for resource type"

  override protected def invalidValueMessage(value: WdlValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object InstanceTypeValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new InstanceTypeValidation(config).optional
}

class InstanceTypeValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("instanceType") with OptionalWithDefault[String] {
  override protected def usedInCallCaching: Boolean = false

  override protected def missingValueMessage: String = "Can't find an attribute value for key instance type"

  override protected def invalidValueMessage(value: WdlValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object ImageIdValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new ImageIdValidation(config).optional
}

class ImageIdValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("imageId") with OptionalWithDefault[String] {
  override protected def usedInCallCaching: Boolean = false

  override protected def missingValueMessage: String = "Can't find an attribute value for key image id"

  override protected def invalidValueMessage(value: WdlValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object WorkerPathValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new WorkerPathValidation(config).optional
}

class WorkerPathValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("workerPath") with OptionalWithDefault[String] {
  override protected def usedInCallCaching: Boolean = false

  override protected def missingValueMessage: String = "Can't find an attribute value for key worker path"

  override protected def invalidValueMessage(value: WdlValue): String = super.missingValueMessage

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[String]] = {
    case WdlString(value) => value.validNel
  }
}

object ReserveOnFailValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = new ReserveOnFailValidation(config).optional
}

class ReserveOnFailValidation(override val config: Option[Config]) extends BooleanRuntimeAttributesValidation(BcsRuntimeAttributes.ReserveOnFailKey) with OptionalWithDefault[Boolean]

object AutoReleaseJobValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = new AutoReleaseJobValidation(config).optional
}

class AutoReleaseJobValidation(override val config: Option[Config]) extends BooleanRuntimeAttributesValidation(BcsRuntimeAttributes.AutoReleaseJobKey) with OptionalWithDefault[Boolean]

object TimeoutValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Int] = new TimeoutValidation(config).optional
}

class TimeoutValidation(override val config: Option[Config]) extends IntRuntimeAttributesValidation(BcsRuntimeAttributes.TimeoutKey) with OptionalWithDefault[Int]

object VerboseValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = new VerboseValidation(config).optional
}

class VerboseValidation(override val config: Option[Config]) extends BooleanRuntimeAttributesValidation(BcsRuntimeAttributes.VerboseKey) with OptionalWithDefault[Boolean]
