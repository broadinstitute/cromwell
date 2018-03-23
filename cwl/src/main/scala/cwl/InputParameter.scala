package cwl

import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import shapeless.Poly1
import wom.callable.Callable.InputDefinition.InputValueMapper
import wom.expression.IoFunctionSet
import wom.types.{WomSingleFileType, WomType}
import wom.values.{WomArray, WomMaybePopulatedFile, WomObject, WomObjectLike, WomOptionalValue, WomValue}

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

  /**
    * Yet another value mapper. This one is needed because in CWL we might need to "augment" inputs which we can only do
    * once they have been linked to a WomValue. This input value mapper encapsulates logic to be applied once that is
    * done. For now, if the inputParameter has an input binding with loadContents = true, load the content of the file.
    *
    * This is based on the spec in http://www.commonwl.org/v1.0/CommandLineTool.html#Input_binding
    *
    * NOTE: There may be _many_ cases not implemented here that need to be fixed.
    */
  def inputValueMapper(inputParameter: InputParameter,
                       inputType: MyriadInputType,
                       expressionLib: ExpressionLib): InputValueMapper = {
    ioFunctionSet: IoFunctionSet => {

      def populateFiles(womValue: WomValue): ErrorOr[WomValue] = {
        womValue match {
          case womMaybePopulatedFile: WomMaybePopulatedFile =>
            val parameterContext = ParameterContext(self = womMaybePopulatedFile)
            val secondaryFilesFromInputParameter = inputParameter.secondaryFiles
            val secondaryFilesFromType = inputType.fold(MyriadInputTypeToSecondaryFiles)
            val secondaryFiles = secondaryFilesFromInputParameter orElse secondaryFilesFromType
            for {
              loaded <- maybeLoadContents(womMaybePopulatedFile, ioFunctionSet, inputParameter.loadContents)
              secondaries <- FileParameter.secondaryFiles(
                loaded,
                WomSingleFileType,
                secondaryFiles,
                parameterContext,
                expressionLib
              )
              updated = loaded.copy(secondaryFiles = secondaries)
            } yield updated

          case WomArray(_, values) => values.toList.traverse(populateFiles).map(WomArray(_))
          case WomOptionalValue(_, Some(innerValue)) => populateFiles(innerValue).map(WomOptionalValue(_))
          case obj: WomObjectLike =>
            // Map the values
            obj.values.toList.traverse[ErrorOr, (String, WomValue)]({
              case (key, value) => populateFiles(value).map(key -> _)
            })
              .map(_.toMap)
              // transform to Either so we can flatMap
              .toEither
              // Validate new types are still valid w.r.t the object
              .flatMap(WomObject.withTypeChecked(_, obj.womObjectTypeLike))
              // re-transform to ErrorOr
              .toValidated
          case womValue: WomValue =>
            womValue.valid
        }
      }

      womValue => populateFiles(womValue).toTry(s"loading $womValue for ${inputParameter.id}").get
    }
  }

  private def maybeLoadContents(womMaybePopulatedFile: WomMaybePopulatedFile,
                                ioFunctionSet: IoFunctionSet,
                                loadContents: Boolean): ErrorOr[WomMaybePopulatedFile] = {
    if (loadContents) {
      FileParameter.load64KiB(womMaybePopulatedFile.value, ioFunctionSet) map { contents =>
        womMaybePopulatedFile.copy(contentsOption = Option(contents))
      }
    } else {
      womMaybePopulatedFile.valid
    }
  }
}
