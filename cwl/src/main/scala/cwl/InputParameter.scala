package cwl

import cats.instances.option._
import cats.syntax.functor._
import cats.syntax.parallel._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.IOChecked._
import common.validation.Validation._
import cwl.FileParameter._
import cwl.ontology.Schema
import shapeless.Poly1
import wom.callable.Callable.InputDefinition.InputValueMapper
import wom.expression.IoFunctionSet
import wom.types.{WomSingleFileType, WomType}
import wom.values.{WomArray, WomMaybeListedDirectory, WomMaybePopulatedFile, WomObject, WomObjectLike, WomOptionalValue, WomValue}

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
              .traverse(_.fold(this).apply(womType))
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

  object InputParameterFormatPoly extends Poly1 {
    implicit val caseExpression: Case.Aux[Expression, ParameterContext => ErrorOr[List[String]]] = {
      at { expression =>
        parameterContext =>
          ExpressionEvaluator.eval(expression, parameterContext) map {
            case WomArray(_, values) => values.toList.map(_.valueString)
            case womValue => List(womValue.valueString)
          }
      }
    }

    implicit val caseString: Case.Aux[String, ParameterContext => ErrorOr[List[String]]] = {
      at { string =>
        _ =>
          List(string).valid
      }
    }

    implicit val caseArrayString: Case.Aux[Array[String], ParameterContext => ErrorOr[List[String]]] = {
      at { array =>
        _ =>
          array.toList.valid
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
                       expressionLib: ExpressionLib,
                       schemaOption: Option[Schema]): InputValueMapper = {
    ioFunctionSet: IoFunctionSet => {
      import ioFunctionSet.cs

      def populateFiles(womValue: WomValue): IOChecked[WomValue] = {
        womValue match {
          case womMaybePopulatedFile: WomMaybePopulatedFile =>
            // Don't include the secondary files in the self variables
            val parameterContext = ParameterContext(ioFunctionSet, expressionLib, self = womMaybePopulatedFile.copy(secondaryFiles = List.empty))
            val secondaryFilesFromInputParameter = inputParameter.secondaryFiles
            val secondaryFilesFromType = inputType.fold(MyriadInputTypeToSecondaryFiles)
            val secondaryFiles = secondaryFilesFromInputParameter orElse secondaryFilesFromType
            val inputFormatsErrorOr = inputParameter.format
                .traverse(_.fold(InputParameterFormatPoly).apply(parameterContext))

            for {
              inputFormatsOption <- inputFormatsErrorOr.toIOChecked
              _ <- checkFormat(womMaybePopulatedFile, inputFormatsOption, schemaOption).toIOChecked
              contentsOption <- FileParameter.maybeLoadContents(
                womMaybePopulatedFile,
                ioFunctionSet,
                inputParameter.loadContents
              ).to[IOChecked]
              withSize <- womMaybePopulatedFile.withSize(ioFunctionSet).to[IOChecked]
              loaded = withSize.copy(contentsOption = contentsOption)
              secondaries <- FileParameter.secondaryFiles(
                loaded,
                WomSingleFileType,
                secondaryFiles,
                parameterContext,
                expressionLib,
                ioFunctionSet
              )
              updated = loaded.copy(secondaryFiles = loaded.secondaryFiles ++  secondaries)
            } yield updated
          case womMaybeListedDirectory: WomMaybeListedDirectory => womMaybeListedDirectory.withSize(ioFunctionSet).to[IOChecked].widen
          case WomArray(_, values) => values.toList.parTraverse[IOChecked, IOCheckedPar, WomValue](populateFiles).map(WomArray(_))
          case WomOptionalValue(_, Some(innerValue)) => populateFiles(innerValue).map(WomOptionalValue(_))
          case obj: WomObjectLike =>
            // Map the values
            val populated: IOChecked[WomObject] = obj.values.toList.parTraverse[IOChecked, IOCheckedPar, (String, WomValue)]({
              case (key, value) => populateFiles(value).map(key -> _)
            })
            .map(_.toMap)
            // Validate new types are still valid w.r.t the object
            .flatMap(WomObject.withTypeChecked(_, obj.womObjectTypeLike).toIOChecked)
            
            populated.widen[WomValue]
          case womValue: WomValue =>
            pure(womValue)
        }
      }

      womValue => populateFiles(womValue)
    }
  }
}
