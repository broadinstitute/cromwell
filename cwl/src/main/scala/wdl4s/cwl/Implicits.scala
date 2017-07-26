package wdl4s.cwl

import io.circe.Decoder
import io.circe._
import io.circe.generic.auto._

import io.circe.shapes._
import io.circe.generic.auto._
import shapeless.Coproduct
import cats.syntax.either._
import eu.timepit.refined.string._
import eu.timepit.refined._
import io.circe.refined._
import io.circe.literal._
import io.circe.syntax._
import cats.data.ValidatedNel
import cats.data.Validated._
import cats.syntax.traverse._
import cats.instances.list._
import cats.syntax.option._
import lenthall.validation.ErrorOr.ErrorOr

object Implicits {

  implicit val cwlTypeDecoder = Decoder.enumDecoder(CwlType)
  implicit val cwlVersionDecoder = Decoder.enumDecoder(CwlVersion)
  implicit val scatterMethodDecoder = Decoder.enumDecoder(ScatterMethod)
  implicit val linkMergeMethodDecoder = Decoder.enumDecoder(LinkMergeMethod)

  implicit val cwlTypeEncoder = Encoder.enumEncoder(CwlType)
  implicit val cwlVersionEncoder = Encoder.enumEncoder(CwlVersion)
  implicit val scatterMethodEncoder = Encoder.enumEncoder(ScatterMethod)
  implicit val linkMergeMethodEncoder = Encoder.enumEncoder(LinkMergeMethod)

  implicit val envvarD = Decoder[EVR]
  implicit val ijr = Decoder[IJR]
  implicit val sr = Decoder[SR]
  implicit val sw = Decoder[SWFR]
  implicit def sdr = Decoder[SDR]
  implicit def rr = Decoder[RR]
  implicit def iwdr = Decoder[IWDR]
  implicit val dr = Decoder[DR]
  implicit val scr = Decoder[SCR]
  implicit val sfr = Decoder[SFR]
  implicit val mifr = Decoder[MIFR]
  implicit val sier = Decoder[SIER]

  def s[T <: Function1[S, U], S <: String, U](t: Target, s: S)(implicit
    selector: shapeless.ops.coproduct.Selector[Target,T],
    inj: shapeless.ops.coproduct.Inject[wdl4s.cwl.Requirement, U]) =
      t.select[T].toValidNel(s"Expecting a $s but got $t instead.").map(_(s)).map(Coproduct[Requirement](_))

  def select[T](t: Target)(implicit selector: shapeless.ops.coproduct.Selector[Target,T]):ValidatedNel[String, T] =
             t.select[T].toValidNel(s"Expecting a EnvVarRequirement but got $t instead.")

  implicit val requirementArrayMapParser: Decoder[Array[Requirement]] = {

    def req[T](a:ErrorOr[T])(implicit inj: shapeless.ops.coproduct.Inject[wdl4s.cwl.Requirement, T]) = a.map(Coproduct[Requirement](_))


    Decoder[Map[String, Target]].
      emap {
        _.toList.traverse[ErrorOr,Requirement] {
            case ("EnvVarRequirement", target) => s[EVR, W.`"EnvVarRequirement"`.T, EnvVarRequirement](target,"EnvVarRequirement")
            //TODO: clean these up to use this approach
            //case ("EnvVarRequirement", target) => TargetFunction[EnvVarRequirement].apply(target)
            case ("InlineJavascriptRequirement", target) => req(select[IJR](target).map(_("InlineJavascriptRequirement")))
            case ("SchemaDefRequirement", target) => req(select[SDR](target).map(_("SchemaDefRequirement")))
            case ("DockerRequirement", target) => req(select[DR](target).map(_("DockerRequirement")))
            case ("SoftwareRequirement", target) => req(select[SR](target).map(_("SoftwareRequirement")))
            case ("InitialWorkDirRequirement", target) => req(select[IWDR](target).map(_("InitialWorkDirRequirement")))
            case ("ShellCommandRequirement", target) => req(select[SCR](target).map(_("ShellCommandRequirement")))
            case ("ResourceRequirement", target) => req(select[RR](target).map(_("ResourceRequirement")))
            case ("SubworkflowFeatureRequirement", target) => req(select[SWFR](target).map(_("SubworkflowFeatureRequirement")))
            case ("ScatterFeatureRequirement", target) => req(select[SFR](target).map(_("ScatterFeatureRequirement")))
            case ("MultipleInputFeatureRequirement", target) => req(select[MIFR](target).map(_("MultipleInputFeatureRequirement")))
            case ("StepInputExpressionRequirement", target) => req(select[SIER](target).map(_("StepInputExpressionRequirement")))
            case (key,_) => invalidNel(s"key $key was not amongst possible values " +
              "InlineJavascriptRequirement, SchemaDefRequirement, DockerRequirement, SoftwareRequirement, InitialWorkDirRequirement, EnvVarRequirement, ShellCommandRequirement, ResourceRequirement, SubworkflowFeatureRequirement, ScatterFeatureRequirement, MultipleInputFeatureRequirement, StepInputExpressionRequirement")
          }.toEither.leftMap(_.toList.mkString(", ")).map(_.toArray)
      }
  }

  implicit def enumerationEncoder[V <: Enumeration#Value]: Encoder[V] = (value: V) => Json.fromString(value.toString)
}
