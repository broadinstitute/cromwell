package cwl

import cats.syntax.validated._
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.Validation._
import common.validation.IOChecked._
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import cwl.InitialWorkDirFileGeneratorExpression._
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import shapeless.Poly1
import wom.callable.{AdHocValue, ContainerizedInputExpression}
import wom.expression.IoFunctionSet.{IoDirectory, IoFile}
import wom.expression.{FileEvaluation, IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

import scala.concurrent.Await
import scala.concurrent.duration.Duration

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

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = Set.empty[FileEvaluation].validNel
}

final case class InitialWorkDirFileGeneratorExpression(entry: IwdrListingArrayEntry, expressionLib: ExpressionLib) extends ContainerizedInputExpression {

  def evaluate(inputValues: Map[String, WomValue], mappedInputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): IOChecked[List[AdHocValue]] = {
    def recursivelyBuildDirectory(directory: String): IOChecked[WomMaybeListedDirectory] = {
      import cats.instances.list._
      import cats.syntax.traverse._
      for {
        listing <- ioFunctionSet.listDirectory(directory)().toIOChecked
        fileListing <- listing.toList.traverse[IOChecked, WomFile] {
          case IoDirectory(e) => recursivelyBuildDirectory(e).widen
          case IoFile(e) => WomSingleFile(e).validIOChecked.widen
        }
      } yield WomMaybeListedDirectory(Option(directory), Option(fileListing))
    }
    
    inputValues.toList.traverse[IOChecked, (String, WomValue)]({
      case (k, v: WomMaybeListedDirectory) =>
        val absolutePathString = ioFunctionSet.pathFunctions.relativeToHostCallRoot(v.value)
        recursivelyBuildDirectory(absolutePathString).contextualizeErrors(s"Error building directory $absolutePathString") map { k -> _ }
      case kv => kv.validIOChecked
    }).map(_.toMap)
      .flatMap({ updatedValues =>
        val unmappedParameterContext = ParameterContext(ioFunctionSet, expressionLib, updatedValues)
        entry.fold(InitialWorkDirFilePoly).apply(unmappedParameterContext, mappedInputValues).toIOChecked
      })
  }
}

object InitialWorkDirFileGeneratorExpression {
  type InitialWorkDirFileEvaluator = (ParameterContext, Map[String, WomValue]) => ErrorOr[List[AdHocValue]]

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

          val womValueErrorOr: ErrorOr[AdHocValue] = entryEvaluation flatMap {
            case womFile: WomFile =>
              val errorOrEntryName: ErrorOr[Option[String]] = expressionDirent.entryname match {
                case Some(actualEntryName) => actualEntryName.fold(EntryNamePoly).apply(unmappedParameterContext).map(Option.apply)
                case None => None.valid
              }
              errorOrEntryName map { entryName =>
                AdHocValue(womFile, entryName, inputName = mutableInputOption)
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
            } yield AdHocValue(writtenFile, alternativeName = None, inputName = mutableInputOption)
          }

          womValueErrorOr.map(List(_))
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

            womValueErrorOr.map(AdHocValue(_, alternativeName = None, inputName = None)).map(List(_))
        }
      }
    }

    implicit val caseExpression: Case.Aux[Expression, InitialWorkDirFileEvaluator] = {
      at { expression =>
        (unmappedParameterContext, _) => {
          // A single expression which must evaluate to an array of Files
          val expressionEvaluation = ExpressionEvaluator.eval(expression, unmappedParameterContext)

          expressionEvaluation flatMap {
            case array: WomArray if array.value.forall(_.isInstanceOf[WomFile]) => 
              array.value.toList.map(_.asInstanceOf[WomFile]).map(AdHocValue(_, alternativeName = None, inputName = None)).validNel
            case file: WomFile =>
              List(AdHocValue(file, alternativeName = None, inputName = None)).validNel
            case other =>
              val error = "InitialWorkDirRequirement listing expression must be File or Array[File] but got %s: %s"
                .format(other, other.womType.stableName)
              error.invalidNel
          }
        }
      }
    }

    implicit val caseString: Case.Aux[String, InitialWorkDirFileEvaluator] = {
      at { string =>
        (_, _) => {
          List(AdHocValue(WomSingleFile(string), alternativeName = None, inputName = None)).valid
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
          file.asWomValue.map(AdHocValue(_, alternativeName = None, inputName = None)).map(List(_))
        }
      }
    }

    implicit val caseDirectory: Case.Aux[Directory, InitialWorkDirFileEvaluator] = {
      at { directory =>
        (_, _) => {
          directory.asWomValue.map(AdHocValue(_, alternativeName = None, inputName = None)).map(List(_))
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
