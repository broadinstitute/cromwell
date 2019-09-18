package cromwell.backend.impl.bcs

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import cromwell.backend.impl.bcs.BcsClusterIdOrConfiguration.BcsClusterIdOrConfiguration
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation._
import net.ceedubs.ficus.Ficus._
import wom.types._
import wom.values._

import scala.util.{Failure, Success, Try}


trait OptionalWithDefault[A] {
  this: RuntimeAttributesValidation[A] =>
  protected val config: Option[Config]

  override protected def staticDefaultOption: Option[WomValue] = {
    Try(this.configDefaultWomValue(config)) match {
      case Success(value: Option[WomValue]) => value
      case Failure(_) => None
    }
  }
}

final case class BcsRuntimeAttributes(continueOnReturnCode: ContinueOnReturnCode,
                                dockerTag: Option[BcsDocker],
                                docker: Option[BcsDocker],
                                failOnStderr: Boolean,
                                mounts: Option[Seq[BcsMount]],
                                userData: Option[Seq[BcsUserData]],
                                cluster: Option[BcsClusterIdOrConfiguration],
                                imageId: Option[String],
                                systemDisk: Option[BcsSystemDisk],
                                dataDisk: Option[BcsDataDisk],
                                reserveOnFail: Option[Boolean],
                                autoReleaseJob: Option[Boolean],
                                timeout: Option[Int],
                                verbose: Option[Boolean],
                                vpc: Option[BcsVpcConfiguration],
                                tag: Option[String])

object BcsRuntimeAttributes {

  val MountsKey = "mounts"
  val UserDataKey = "userData"
  val MountsDefaultValue = WomString("")
  val ReserveOnFailKey = "reserveOnFail"
  val ReserveOnFailDefault = false
  val AutoReleaseJobKey = "autoReleaseJob"
  val AutoReleaseJobDefault = WomBoolean(true)
  val TimeoutKey = "timeout"
  val TimeoutDefault = WomInteger(21600)
  val VerboseKey = "verbose"
  val ClusterKey = "cluster"
  val DockerKey = "docker"
  val SystemDiskKey = "systemDisk"
  val DataDiskKey = "dataDisk"
  val VpcKey = "vpc"
  val TagKey = "tag"

  private def failOnStderrValidation(runtimeConfig: Option[Config]) = FailOnStderrValidation.default(runtimeConfig)

  private def continueOnReturnCodeValidation(runtimeConfig: Option[Config]) = ContinueOnReturnCodeValidation.default(runtimeConfig)

  private def clusterValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsClusterIdOrConfiguration] = ClusterValidation.optionalWithDefault(runtimeConfig)

  private def dockerTagValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsDocker] = DockerTagValidation.optionalWithDefault(runtimeConfig)
  private def dockerValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsDocker] = DockerValidation.optionalWithDefault(runtimeConfig)

  private def userDataValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsUserData]] = UserDataValidation.optionalWithDefault(runtimeConfig)

  private def systemDiskValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsSystemDisk] = SystemDiskValidation.optionalWithDefault(runtimeConfig)
  private def dataDiskValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsDataDisk] = DataDiskValidation.optionalWithDefault(runtimeConfig)

  private def reserveOnFailValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = ReserveOnFailValidation.optionalWithDefault(runtimeConfig)

  private def autoReleaseJobValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = AutoReleaseJobValidation.optionalWithDefault(runtimeConfig)

  private def mountsValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsMount]] = MountsValidation.optionalWithDefault(runtimeConfig)

  private def timeoutValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Int] = TimeoutValidation.optionalWithDefault(runtimeConfig)

  private def verboseValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[Boolean] = VerboseValidation.optionalWithDefault(runtimeConfig)

  private def vpcValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[BcsVpcConfiguration] = VpcValidation.optionalWithDefault(runtimeConfig)

  private def tagValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = TagValidation.optionalWithDefault(runtimeConfig)

  private def imageIdValidation(runtimeConfig: Option[Config]): OptionalRuntimeAttributesValidation[String] = ImageIdValidation.optionalWithDefault(runtimeConfig)

  def runtimeAttributesBuilder(backendRuntimeConfig: Option[Config]): StandardValidatedRuntimeAttributesBuilder = {
    val defaults = StandardValidatedRuntimeAttributesBuilder.default(backendRuntimeConfig).withValidation(
      mountsValidation(backendRuntimeConfig),
      userDataValidation(backendRuntimeConfig),
      clusterValidation(backendRuntimeConfig),
      systemDiskValidation(backendRuntimeConfig),
      dataDiskValidation(backendRuntimeConfig),
      reserveOnFailValidation(backendRuntimeConfig),
      autoReleaseJobValidation(backendRuntimeConfig),
      timeoutValidation(backendRuntimeConfig),
      verboseValidation(backendRuntimeConfig),
      vpcValidation(backendRuntimeConfig),
      tagValidation(backendRuntimeConfig),
      imageIdValidation(backendRuntimeConfig)
    )

    // TODO: docker trips up centaur testing, for now https://github.com/broadinstitute/cromwell/issues/3518
    if (backendRuntimeConfig.exists(_.getOrElse("ignoreDocker", false))) {
      defaults
    } else {
      defaults.withValidation(
        dockerTagValidation(backendRuntimeConfig),
        dockerValidation(backendRuntimeConfig)
      )
    }
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, backendRuntimeConfig: Option[Config]): BcsRuntimeAttributes = {
    val failOnStderr: Boolean =
      RuntimeAttributesValidation.extract(failOnStderrValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val continueOnReturnCode: ContinueOnReturnCode =
      RuntimeAttributesValidation.extract(continueOnReturnCodeValidation(backendRuntimeConfig), validatedRuntimeAttributes)
    val mounts: Option[Seq[BcsMount]] = RuntimeAttributesValidation.extractOption(mountsValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val userData: Option[Seq[BcsUserData]] = RuntimeAttributesValidation.extractOption(userDataValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    val cluster: Option[BcsClusterIdOrConfiguration] = RuntimeAttributesValidation.extractOption(clusterValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val imageId: Option[String] = RuntimeAttributesValidation.extractOption(imageIdValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val dockerTag: Option[BcsDocker] = RuntimeAttributesValidation.extractOption(dockerTagValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val docker: Option[BcsDocker] = RuntimeAttributesValidation.extractOption(dockerValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val systemDisk: Option[BcsSystemDisk] = RuntimeAttributesValidation.extractOption(systemDiskValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val dataDisk: Option[BcsDataDisk] = RuntimeAttributesValidation.extractOption(dataDiskValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    val reserveOnFail: Option[Boolean] = RuntimeAttributesValidation.extractOption(reserveOnFailValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val autoReleaseJob: Option[Boolean] = RuntimeAttributesValidation.extractOption(autoReleaseJobValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val timeout: Option[Int] = RuntimeAttributesValidation.extractOption(timeoutValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val verbose: Option[Boolean] = RuntimeAttributesValidation.extractOption(verboseValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val vpc: Option[BcsVpcConfiguration] = RuntimeAttributesValidation.extractOption(vpcValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)
    val tag: Option[String] = RuntimeAttributesValidation.extractOption(tagValidation(backendRuntimeConfig).key, validatedRuntimeAttributes)

    new BcsRuntimeAttributes(
      continueOnReturnCode,
      dockerTag,
      docker,
      failOnStderr,
      mounts,
      userData,
      cluster,
      imageId,
      systemDisk,
      dataDisk,
      reserveOnFail,
      autoReleaseJob,
      timeout,
      verbose,
      vpc,
      tag
    )
  }
}

object MountsValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsMount]] = new MountsValidation(config).optional
}

class MountsValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[Seq[BcsMount]] with OptionalWithDefault[Seq[BcsMount]] {
  override def key: String = BcsRuntimeAttributes.MountsKey

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[BcsMount]]] = {
    case WomString(value) => validateMounts(value.split(",\\s*").toSeq)
    case WomArray(wdlType, values) if wdlType.memberType == WomStringType =>
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
    val emptyMountNel: ErrorOr[Vector[BcsMount]] = Vector.empty[BcsMount].validNel
    val mountsNel: ErrorOr[Vector[BcsMount]] = nels.foldLeft(emptyMountNel) {
      (acc, v) => (acc, v) mapN { (a, v) => a :+ v }
    }
    mountsNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}

object UserDataValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[Seq[BcsUserData]] = new UserDataValidation(config).optional
}

class UserDataValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[Seq[BcsUserData]] with OptionalWithDefault[Seq[BcsUserData]]{
  override def key: String = BcsRuntimeAttributes.UserDataKey

  override def usedInCallCaching: Boolean = true

  override def coercion: Traversable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[BcsUserData]]] = {
    case WomString(value) => validateUserData(value.split(",\\s*").toSeq)
    case WomArray(wdlType, values) if wdlType.memberType == WomStringType =>
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
    val emptyDataNel: ErrorOr[Vector[BcsUserData]] = Vector.empty[BcsUserData].validNel
    val datasNel: ErrorOr[Vector[BcsUserData]] = nels.foldLeft(emptyDataNel) {
      (acc, v) => (acc, v) mapN { (a, v) => a :+ v }
    }
    datasNel
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
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

  override def coercion: Traversable[WomType] = Set(WomStringType)

  override def validateValue: PartialFunction[WomValue, ErrorOr[BcsClusterIdOrConfiguration]] = {
    case WomString(s) => BcsClusterIdOrConfiguration.parse(s.toString) match {
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
  override def coercion: Traversable[WomType] = Set(WomStringType)
  override def validateValue: PartialFunction[WomValue, ErrorOr[BcsSystemDisk]] = {
    case WomString(s) => BcsDisk.parse(s.toString) match {
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
  override def coercion: Traversable[WomType] = Set(WomStringType)
  override def validateValue: PartialFunction[WomValue, ErrorOr[BcsDataDisk]] = {
    case WomString(s) => BcsDisk.parse(s.toString) match {
      case Success(disk: BcsDataDisk) => disk.validNel
      case _ => s"system disk should be string like 'cloud 40 /home/data/'".invalidNel
    }
  }
}

object DockerTagValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsDocker] = new DockerTagValidation(config).optional
}

class DockerTagValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[BcsDocker] with OptionalWithDefault[BcsDocker]
{
  override def key: String = "dockerTag"
  override def coercion: Traversable[WomType] = Set(WomStringType)
  override def validateValue: PartialFunction[WomValue, ErrorOr[BcsDocker]] = {
    case WomString(s) => BcsDocker.parse(s.toString) match {
      case Success(docker: BcsDocker) => docker.validNel
      case _ => s"docker must be 'dockerImage dockerPath' like".invalidNel
    }
  }
}

object DockerValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsDocker] = new DockerValidation(config).optional
}

class DockerValidation(override val config: Option[Config]) extends DockerTagValidation(config)
{
  override def key: String = "docker"
  override def usedInCallCaching: Boolean = true
}

object VpcValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[BcsVpcConfiguration] = new VpcValidation(config).optional
}

class VpcValidation(override val config: Option[Config]) extends RuntimeAttributesValidation[BcsVpcConfiguration] with OptionalWithDefault[BcsVpcConfiguration]
{
  override def key: String = "vpc"
  override def coercion: Traversable[WomType] = Set(WomStringType)
  override def validateValue: PartialFunction[WomValue, ErrorOr[BcsVpcConfiguration]] = {
    case WomString(s) => BcsVpcConfiguration.parse(s.toString) match {
      case Success(vpc: BcsVpcConfiguration) => vpc.validNel
      case _ => s"vpc must be '192.168.0.0/16 vpc-xxxx' like".invalidNel
    }
  }
}

object TagValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new TagValidation(config).optional
}

class TagValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("tag") with OptionalWithDefault[String]

object ImageIdValidation {
  def optionalWithDefault(config: Option[Config]): OptionalRuntimeAttributesValidation[String] = new ImageIdValidation(config).optional
}

class ImageIdValidation(override val config: Option[Config]) extends StringRuntimeAttributesValidation("imageId") with OptionalWithDefault[String]
{
  override def usedInCallCaching: Boolean = true
}

