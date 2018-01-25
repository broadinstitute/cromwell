package cwl

import cats.syntax.traverse._
import common.validation.ErrorOr._
import common.validation.Validation._
import shapeless.Poly1
import wom.types.WomType
import wom.values.{WomArray, WomValue}

trait InputParameter {
  def id: String
  def label: Option[String]
  def secondaryFiles: Option[SecondaryFiles]
  def format: Option[InputParameterFormat]
  def streamable: Option[Boolean]
  def doc: Option[Doc]
  def inputBinding: Option[InputCommandLineBinding]
  def default: Option[CwlAny]
  def `type`: Option[MyriadInputType]

  def loadContents = inputBinding.flatMap(_.loadContents).getOrElse(false)
}

object InputParameter {
  object IdDefaultAndType {
    def unapply(arg: InputParameter): Option[(String, CwlAny, MyriadInputType)] = (arg.default, arg.`type`) match {
      case (Some(default), Some(tpe)) => Option((arg.id, default, tpe))
      case _ => None
    }
  }

  object IdAndType {
    def unapply(arg: InputParameter): Option[(String, MyriadInputType)] = (arg.default, arg.`type`) match {
      case (None, Some(tpe)) => Option((arg.id, tpe))
      case _ => None
    }
  }

  type DefaultToWomValueFunction = WomType => ErrorOr[WomValue]

  import cats.instances.list._

  object DefaultToWomValuePoly extends Poly1 {
    implicit def caseFileOrDirectory: Case.Aux[FileOrDirectory, DefaultToWomValueFunction] = {
      at {
        _.fold(this)
      }
    }

    implicit def caseFileOrDirectoryArray: Case.Aux[Array[FileOrDirectory], DefaultToWomValueFunction] = {
      at {
        fileOrDirectoryArray =>
          womType =>
            fileOrDirectoryArray
              .toList
              .traverse[ErrorOr, WomValue](_.fold(this).apply(womType))
              .map(WomArray(_))
      }
    }

    implicit def caseFile: Case.Aux[File, DefaultToWomValueFunction] = {
      at {
        file =>
          womType =>
            file.asWomValue.flatMap(womType.coerceRawValue(_).toErrorOr)
      }
    }

    implicit def caseDirectory: Case.Aux[Directory, DefaultToWomValueFunction] = {
      at {
        directory =>
          womType =>
            directory.asWomValue.flatMap(womType.coerceRawValue(_).toErrorOr)
      }
    }

    implicit def caseJson: Case.Aux[io.circe.Json, DefaultToWomValueFunction] = {
      at {
        circeJson =>
          womType =>
            val stringJson = circeJson.noSpaces
            import spray.json._
            val sprayJson = stringJson.parseJson
            womType.coerceRawValue(sprayJson).toErrorOr
      }
    }
  }
}
