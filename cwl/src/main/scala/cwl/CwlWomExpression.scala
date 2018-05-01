package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.Validation._
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import cwl.FileParameter.sync
import cwl.InitialWorkDirFileGeneratorExpression._
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import shapeless.Poly1
import wom.callable.{AdHocValue, ContainerizedInputExpression}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, Future}

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

  def evaluate(inputValues: Map[String, WomValue], mappedInputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[AdHocValue] = {
    def recursivelyBuildDirectory(directory: String): ErrorOr[WomMaybeListedDirectory] = {
      import cats.instances.list._
      import cats.syntax.traverse._
      for {
        listing <- sync(ioFunctionSet.listDirectory(directory)).toErrorOr
        fileListing <- listing.toList.traverse[ErrorOr, WomFile] {
          case e if sync(ioFunctionSet.isDirectory(e)).getOrElse(false) => recursivelyBuildDirectory(e)
          case e => WomSingleFile(e).validNel
        }
      } yield WomMaybeListedDirectory(Option(directory), Option(fileListing))
    }
    val updatedValues = inputValues map {
      case (k, v: WomMaybeListedDirectory) => k -> recursivelyBuildDirectory(v.value).getOrElse(throw new RuntimeException("boom"))
      case (k, v: WomMaybePopulatedFile) => k -> WomSingleFile(v.value)
      case kv => kv
    }
    val unmappedParameterContext = ParameterContext(ioFunctionSet, expressionLib, updatedValues)
    entry.fold(InitialWorkDirFilePoly).apply(unmappedParameterContext, mappedInputValues)
  }
}

object InitialWorkDirFileGeneratorExpression {
  type InitialWorkDirFileEvaluator = (ParameterContext, Map[String, WomValue]) => ErrorOr[AdHocValue]

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

          val mutableInputOption = expressionDirent.entry.fold(ExpressionToMutableInputOptionPoly)

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

          val womValueErrorOr: ErrorOr[WomSingleFile] = entryEvaluation flatMap {
            case womFile: WomFile =>
              val errorOrEntryName: ErrorOr[String] = expressionDirent.entryname match {
                case Some(actualEntryName) => actualEntryName.fold(EntryNamePoly).apply(unmappedParameterContext)
                case None => unmappedParameterContext.ioFunctionSet.pathFunctions.name(womFile.value).valid
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

          womValueErrorOr.map(AdHocValue(_, mutableInputOption))
        }
      }
    }

    implicit val caseStringDirent: Case.Aux[StringDirent, InitialWorkDirFileEvaluator] = {
      at {
        stringDirent => {
          (unmappedParameterContext, _) =>
            val womValueErrorOr = for {
              entryName <- stringDirent.entryname.fold(EntryNamePoly).apply(unmappedParameterContext)
              contentString = stringDirent.entry
              writeFile = unmappedParameterContext.ioFunctionSet.writeFile(entryName, contentString)
              writtenFile <- validate(Await.result(writeFile, Duration.Inf))
            } yield writtenFile

            womValueErrorOr.map(AdHocValue(_, mutableInputOption = None))
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
            val baseFileName = unmappedParameterContext.ioFunctionSet.pathFunctions.name(file.value)
            unmappedParameterContext.ioFunctionSet.copyFile(file.value, baseFileName)
          }

          def stageFiles(womArray: WomArray): ErrorOr[WomValue] = {
            implicit val ec = unmappedParameterContext.ioFunctionSet.ec

            val unstagedFiles = womArray.value.map(_.asInstanceOf[WomFile])
            val stagedFiles = Await.result(Future.sequence(unstagedFiles map stageFile), Duration.Inf)
            womArray.womType.coerceRawValue(stagedFiles).toErrorOr
          }

          val womValueErrorOr = expressionEvaluation flatMap {
            case array: WomArray if array.value.forall(_.isInstanceOf[WomFile]) => stageFiles(array)
            case file: WomFile =>
              validate(Await.result(stageFile(file), Duration.Inf))
            case other =>
              val error = "InitialWorkDirRequirement listing expression must be File or Array[File] but got %s: %s"
                .format(other, other.womType.toDisplayString)
              error.invalidNel
          }

          womValueErrorOr.map(AdHocValue(_, mutableInputOption = None))
        }
      }
    }

    implicit val caseString: Case.Aux[String, InitialWorkDirFileEvaluator] = {
      at { string =>
        (_, _) => {
          AdHocValue(WomSingleFile(string), mutableInputOption = None).valid
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
          file.asWomValue.map(AdHocValue(_, mutableInputOption = None))
        }
      }
    }

    implicit val caseDirectory: Case.Aux[Directory, InitialWorkDirFileEvaluator] = {
      at { directory =>
        (_, _) => {
          directory.asWomValue.map(AdHocValue(_, mutableInputOption = None))
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

  /**
    * Searches for ECMAScript expressions that mutate input paths.
    *
    * @see [[AdHocValue]]
    */
  object ExpressionToMutableInputOptionPoly extends Poly1 {
    private val MutableInputRegex = """(?s)\s*\$\(\s*inputs\.(\w+)\s*\)\s*""".r
    implicit val caseECMAScriptExpression: Case.Aux[ECMAScriptExpression, Option[String]] = {
      at { eCMAScriptExpression =>
        eCMAScriptExpression.value match {
          case MutableInputRegex(mutableInput) => Option(mutableInput)
          case _ => None
        }
      }
    }
    implicit val caseECMAScriptFunction: Case.Aux[ECMAScriptFunction, Option[String]] = at { _ => None }
  }
}
