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
import BcsClusterIdOrConfiguration.BcsClusterIdOrConfiguration

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
                                docker: Option[BcsDocker],
                                failOnStderr: Boolean,
                                mounts: Option[Seq[BcsMount]],
                                userData: Option[Seq[BcsUserData]],
                                cluster: Option[BcsClusterIdOrConfiguration],
                                systemDisk: Option[BcsSystemDisk],
                                dataDisk: Option[BcsDataDisk],
                                reserveOnFail: Option[Boolean],
                                autoReleaseJob: Option[Boolean],
                                workerPath: Option[String],
                                timeout: Option[Int],
                                verbose: Option[Boolean],
                                vpc: Option[BcsVpcConfiguration])

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
  val ClusterKey = "cluster"
  val DockerKey = "docker"
  val SystemDiskKey = "systemDisk"
  val DataDiskKey = "dataDisk"
  val VpcKey = "vpc"

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def clusterValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsClusterIdOrConfiguration] = ClusterValidation.optionalWithDefault(runtimeConfig)

  private def dockerValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsDocker] = DockerValidation.optionalWithDefault(runtimeConfig)
  private def userDataValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsUserData]] = UserDataValidation.optionalWithDefault(runtimeConfig)

  private def systemDiskValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsSystemDisk] = SystemDiskValidation.optionalWithDefault(runtimeConfig)
  private def dataDiskValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsDataDisk] = DataDiskValidation.optionalWithDefault(runtimeConfig)

  private def reserveOnFailValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = ReserveOnFailValidation.optionalWithDefault(runtimeConfig)

  private def autoReleaseJobValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = AutoReleaseJobValidation.optionalWithDefault(runtimeConfig)

  private def mountsValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsMount]] = MountsValidation.optionalWithDefault(runtimeConfig)

  private def workerPathValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = WorkerPathValidation.optionalWithDefault(runtimeConfig)

  private def timeoutValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int] = TimeoutValidation.optionalWithDefault(runtimeConfig)

  private def verboseValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = VerboseValidation.optionalWithDefault(runtimeConfig)

  private def vpcValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsVpcConfiguration] = VpcValidation.optionalWithDefault(runtimeConfig)

  def runtimeAttributesBuilder(backendRuntimeConfig: Option[Config]): StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default(backendRuntimeConfig).withValidation(
      mountsValidation(backendRuntimeConfig),
      userDataValidation(backendRuntimeConfig),
      clusterValidation(backendRuntimeConfig),
      dockerValidation(backendRuntimeConfig),
      systemDiskValidation(backendRuntimeConfig),
      dataDiskValidation(backendRuntimeConfig),
      reserveOnFailValidation(backendRuntimeConfig),
      autoReleaseJobValidation(backendRuntimeConfig),
      workerPathValidation(backendRuntimeConfig),
      timeoutValidation(backendRuntimeConfig),
      verboseValidation(backendRuntimeConfig),
      vpcValidation(backendRuntimeConfig)
    )

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, backendRuntimeConfig: Option[Config]): BcsRuntimeAttributes = {
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val mounts: Option[Seq[BcsMount]] = RuntimeAttributesValidation.extractOption(mountsValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val userData: Option[Seq[BcsUserData]] = RuntimeAttributesValidation.extractOption(userDataValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    val cluster: Option[BcsClusterIdOrConfiguration] = RuntimeAttributesValidation.extractOption(clusterValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val docker: Option[BcsDocker] = RuntimeAttributesValidation.extractOption(dockerValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val systemDisk: Option[BcsSystemDisk] = RuntimeAttributesValidation.extractOption(systemDiskValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val dataDisk: Option[BcsDataDisk] = RuntimeAttributesValidation.extractOption(dataDiskValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    val reserveOnFail: Option[Boolean] = RuntimeAttributesValidation.extractOption(reserveOnFailValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val autoReleaseJob: Option[Boolean] = RuntimeAttributesValidation.extractOption(autoReleaseJobValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val workerPath: Option[String] = RuntimeAttributesValidation.extractOption(workerPathValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val timeout: Option[Int] = RuntimeAttributesValidation.extractOption(timeoutValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val verbose: Option[Boolean] = RuntimeAttributesValidation.extractOption(verboseValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val vpc: Option[BcsVpcConfiguration] = RuntimeAttributesValidation.extractOption(vpcValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    new BcsRuntimeAttributes(
      continueOnReturnCode,
      docker,
      failOnStderr,
      mounts,
      userData,
      cluster,
      systemDisk,
      dataDisk,
      reserveOnFail,
      autoReleaseJob,
      workerPath,
      timeout,
      verbose,
      vpc
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


object ClusterValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsClusterIdOrConfiguration] = new ClusterValidation(config).optional
}

class ClusterValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[BcsClusterIdOrConfiguration] with OptionalWithDefault[BcsClusterIdOrConfiguration]
{
  override def key: String = "cluster"

  override def coercion: Traversable[WdlType] = Set(WdlStringType)

  override def validateValue: PartialFunction[WdlValue, ErrorOr[BcsClusterIdOrConfiguration]] = {
    case WdlString(s) => BcsClusterIdOrConfiguration.parse(s.toString) match {
      case Success(cluster) => cluster.validNel
      case Failure(t) => t.getMessage.invalidNel
    }
  }
}

object SystemDiskValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsSystemDisk] = new SystemDiskValidation(config).optional
}

class SystemDiskValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[BcsSystemDisk] with OptionalWithDefault[BcsSystemDisk]
{
  override def key: String = "systemDisk"
  override def coercion: Traversable[WdlType] = Set(WdlStringType)
  override def validateValue: PartialFunction[WdlValue, ErrorOr[BcsSystemDisk]] = {
    case WdlString(s) => BcsDisk.parse(s.toString) match {
      case Success(disk: BcsSystemDisk) => disk.validNel
      case _ => s"system disk should be string like 'cloud 40'".invalidNel
    }
  }
}

object DataDiskValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsDataDisk] = new DataDiskValidation(config).optional
}

class DataDiskValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[BcsDataDisk] with OptionalWithDefault[BcsDataDisk]
{
  override def key: String = "dataDisk"
  override def coercion: Traversable[WdlType] = Set(WdlStringType)
  override def validateValue: PartialFunction[WdlValue, ErrorOr[BcsDataDisk]] = {
    case WdlString(s) => BcsDisk.parse(s.toString) match {
      case Success(disk: BcsDataDisk) => disk.validNel
      case _ => s"system disk should be string like 'cloud 40 /home/data/'".invalidNel
    }
  }
}

object DockerValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsDocker] = new DockerValidation(config).optional
}

class DockerValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[BcsDocker] with OptionalWithDefault[BcsDocker]
{
  override def key: String = "docker"
  override def coercion: Traversable[WdlType] = Set(WdlStringType)
  override def validateValue: PartialFunction[WdlValue, ErrorOr[BcsDocker]] = {
    case WdlString(s) => BcsDocker.parse(s.toString) match {
      case Success(docker: BcsDocker) => docker.validNel
      case _ => s"docker must be 'dockeImage dockerPath' like".invalidNel
    }
  }
}

object VpcValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsVpcConfiguration] = new VpcValidation(config).optional
}

class VpcValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[BcsVpcConfiguration] with OptionalWithDefault[BcsVpcConfiguration]
{
  override def key: String = "vpc"
  override def coercion: Traversable[WdlType] = Set(WdlStringType)
  override def validateValue: PartialFunction[WdlValue, ErrorOr[BcsVpcConfiguration]] = {
    case WdlString(s) => BcsVpcConfiguration.parse(s.toString) match {
      case Success(vpc: BcsVpcConfiguration) => vpc.validNel
      case _ => s"vpc must be '192.168.0.0/16 vpc-xxxx' like".invalidNel
    }
  }
}