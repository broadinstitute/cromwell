package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.Validation._
import cwl.InitialWorkDirFileGeneratorExpression._
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import shapeless.Poly1
import wom.callable.ContainerizedInputExpression
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, Future}
import scala.util.Try

trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WomType

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = cwlExpressionType.validNel

  def expressionLib: ExpressionLib

  def evaluate(inputs: Map[String, WomValue], parameterContext: ParameterContext, expression: Expression): ErrorOr[WomValue] =
    expression.
      fold(EvaluateExpression).
      apply(parameterContext)
}

case class ECMAScriptWomExpression(expression: Expression,
                                   override val inputs: Set[String],
                                   override val expressionLib: ExpressionLib) extends CwlWomExpression {
  val cwlExpressionType = WomAnyType

  override def sourceString = expression match {
    case Expression.ECMAScriptExpression(s) => s.value
    case Expression.ECMAScriptFunction(s) => s.value
  }

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    val pc = ParameterContext(ioFunctionSet, expressionLib, inputValues)
    evaluate(inputValues, pc, expression)
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = Set.empty[WomFile].validNel
}

final case class InitialWorkDirFileGeneratorExpression(entry: IwdrListingArrayEntry, expressionLib: ExpressionLib) extends ContainerizedInputExpression {

  def evaluate(inputValues: Map[String, WomValue], mappedInputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    val unmappedParameterContext = ParameterContext(ioFunctionSet, expressionLib, inputValues)
    entry.fold(InitialWorkDirFilePoly).apply(unmappedParameterContext, mappedInputValues)
  }
}

object InitialWorkDirFileGeneratorExpression {
  type InitialWorkDirFileEvaluator = (ParameterContext, Map[String, WomValue]) => ErrorOr[WomValue]

  /**
    * Converts an InitialWorkDir.
    *
    * TODO: Review against the spec. Especially for Dirent values. For example:
    *
    * "If the value is an expression that evaluates to a Dirent object, this indicates that the File or Directory in
    * entry should be added to the designated output directory with the name in entryname."
    *
    * - http://www.commonwl.org/v1.0/CommandLineTool.html#InitialWorkDirRequirement
    * - http://www.commonwl.org/v1.0/CommandLineTool.html#Dirent
    */
  object InitialWorkDirFilePoly extends Poly1 {
    implicit val caseExpressionDirent: Case.Aux[ExpressionDirent, InitialWorkDirFileEvaluator] = {
      at { expressionDirent =>
        (unmappedParameterContext, mappedInputValues) => {

          val entryEvaluation =
            expressionDirent.entry match {
                //we need to catch this special case to feed in "value-mapped" input values
              case expr@Expression.ECMAScriptExpression(exprString) if exprString.value.trim() == "$(JSON.stringify(inputs))" =>
                val specialParameterContext = ParameterContext(
                  unmappedParameterContext.ioFunctionSet,
                  unmappedParameterContext.expressionLib,
                  mappedInputValues
                )
                ExpressionEvaluator.eval(expr, specialParameterContext)
              case _ => ExpressionEvaluator.eval(expressionDirent.entry, unmappedParameterContext)
            }

          entryEvaluation flatMap {
            case womFile: WomFile =>
              val errorOrEntryName: ErrorOr[String] = expressionDirent.entryname match {
                case Some(actualEntryName) => actualEntryName.fold(EntryNamePoly).apply(unmappedParameterContext)
                case None => womFile.value.split('/').last.valid
              }
              errorOrEntryName flatMap { entryName =>
                validate {
                  Await.result(unmappedParameterContext.ioFunctionSet.copyFile(womFile.value, entryName), Duration.Inf)
                }
              }
            case other => for {
              coerced <- WomStringType.coerceRawValue(other).toErrorOr
              contentString = coerced.asInstanceOf[WomString].value
              // We force the entryname to be specified, and then evaluate it:
              entryNameStringOrExpression <- expressionDirent.entryname.toErrorOr(
                "Invalid dirent: Entry was a string but no file name was supplied")
              entryName <- entryNameStringOrExpression.fold(EntryNamePoly).apply(unmappedParameterContext)
              writeFile = unmappedParameterContext.ioFunctionSet.writeFile(entryName, contentString)
              writtenFile <- validate(Await.result(writeFile, Duration.Inf))
            } yield writtenFile
          }
        }
      }
    }

    implicit val caseStringDirent: Case.Aux[StringDirent, InitialWorkDirFileEvaluator] = {
      at {
        stringDirent => {
          (unmappedParameterContext, _) =>
            for {
              entryName <- stringDirent.entryname.fold(EntryNamePoly).apply(unmappedParameterContext)
              contentString = stringDirent.entry
              writeFile = unmappedParameterContext.ioFunctionSet.writeFile(entryName, contentString)
              writtenFile <- validate(Await.result(writeFile, Duration.Inf))
            } yield writtenFile
        }
      }
    }

    implicit val caseExpression: Case.Aux[Expression, InitialWorkDirFileEvaluator] = {
      at { expression =>
        (unmappedParameterContext, _) => {
          // A single expression which must evaluate to an array of Files
          val expressionEvaluation = ExpressionEvaluator.eval(expression, unmappedParameterContext)

          def stageFile(file: WomFile): Future[WomSingleFile] = {
            // TODO WomFile could be a WomMaybePopulatedFile with secondary files but this code only stages in
            // the primary file.
            // The file should be staged to the initial work dir using the base filename.
            val baseFileName = file.value.substring(file.value.lastIndexOf('/') + 1)
            unmappedParameterContext.ioFunctionSet.copyFile(file.value, baseFileName)
          }

          def stageFiles(unstagedInputValue: Try[WomValue]): ErrorOr[WomValue] = {
            implicit val ec = unmappedParameterContext.ioFunctionSet.ec

            import common.validation.ErrorOr._
            for {
              womArray <- unstagedInputValue.toErrorOr.map(_.asInstanceOf[WomArray])
              unstagedFiles = womArray.value.map(_.asInstanceOf[WomFile])
              stagedFiles = Await.result(Future.sequence(unstagedFiles map stageFile), Duration.Inf)
              arrayOfStagedFiles <- womArray.womType.coerceRawValue(stagedFiles).toErrorOr
            } yield arrayOfStagedFiles
          }

          import mouse.all._
          expressionEvaluation flatMap {
            case array: WomArray if WomArrayType(WomSingleFileType).coercionDefined(array) =>
              WomArrayType(WomSingleFileType).coerceRawValue(array) |> stageFiles
            case array: WomArray if WomArrayType(WomMaybePopulatedFileType).coercionDefined(array) =>
              WomArrayType(WomMaybePopulatedFileType).coerceRawValue(array) |> stageFiles
            case file: WomFile =>
              validate(Await.result(stageFile(file), Duration.Inf))
            case other =>
              val error = "InitialWorkDirRequirement listing expression must be File or Array[File] but got %s: %s"
                .format(other, other.womType.toDisplayString)
              error.invalidNel
          }
        }
      }
    }

    implicit val caseString: Case.Aux[String, InitialWorkDirFileEvaluator] = {
      at { string =>
        (_, _) => {
          WomSingleFile(string).valid
        }
      }
    }

    implicit val caseStringOrExpression: Case.Aux[StringOrExpression, InitialWorkDirFileEvaluator] = {
      at {
        _.fold(this)
      }
    }

    implicit val caseFile: Case.Aux[File, InitialWorkDirFileEvaluator] = {
      at { file =>
        (_, _) => {
          file.asWomValue
        }
      }
    }

    implicit val caseDirectory: Case.Aux[Directory, InitialWorkDirFileEvaluator] = {
      at { directory =>
        (_, _) => {
          directory.asWomValue
        }
      }
    }

  }

  type EntryNameEvaluator = ParameterContext => ErrorOr[String]

  object EntryNamePoly extends Poly1 {
    implicit val caseString: Case.Aux[String, EntryNameEvaluator] = {
      at {
        string => {
          _ =>
            string.valid
        }
      }
    }

    implicit val caseExpression: Case.Aux[Expression, EntryNameEvaluator] = {
      at {
        expression => {
          parameterContext =>
            for {
              entryNameExpressionEvaluated <- ExpressionEvaluator.eval(expression, parameterContext)
              entryNameValidated <- mustBeString(entryNameExpressionEvaluated)
            } yield entryNameValidated
        }
      }
    }

    private def mustBeString(womValue: WomValue): ErrorOr[String] = {
      womValue match {
        case WomString(s) => s.valid
        case other => WomStringType.coerceRawValue(other).map(_.asInstanceOf[WomString].value).toErrorOr
      }
    }
  }

}
