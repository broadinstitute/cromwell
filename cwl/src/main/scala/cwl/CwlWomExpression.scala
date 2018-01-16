package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.Validation._
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

import scala.concurrent.Await
import scala.concurrent.duration.Duration

trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WomType

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = cwlExpressionType.validNel

  def expressionLib: ExpressionLib

  def evaluate(inputs: Map[String, WomValue], parameterContext: ParameterContext, expression: Expression, expressionLib: ExpressionLib): ErrorOr[WomValue] =
    expression.
      fold(EvaluateExpression).
      apply(parameterContext, expressionLib)
}

case class JobPreparationExpression(expression: Expression,
                                    override val inputs: Set[String],
                                    override val expressionLib: ExpressionLib) extends CwlWomExpression {
  val cwlExpressionType = WomAnyType

  override def sourceString = expression match {
    case Expression.ECMAScriptExpression(s) => s.value
    case Expression.ECMAScriptFunction(s) => s.value
  }

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    val pc = ParameterContext(inputValues)
    evaluate(inputValues, pc, expression, expressionLib)
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = Set.empty[WomFile].validNel
}

final case class InitialWorkDirFileGeneratorExpression(entry: IwdrListingArrayEntry, expressionLib: ExpressionLib) extends CwlWomExpression {
  override def cwlExpressionType: WomType = WomSingleFileType
  override def sourceString: String = entry.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    def mustBeString(womValue: WomValue): ErrorOr[String] = womValue match {
      case WomString(s) => s.validNel
      case other => WomStringType.coerceRawValue(other).map(_.asInstanceOf[WomString].value).toErrorOr
    }

    def evaluateEntryName(stringOrExpression: StringOrExpression): ErrorOr[String] = stringOrExpression match {
      case StringOrExpression.String(s) => s.validNel
      case StringOrExpression.Expression(entrynameExpression) => for {
        entryNameExpressionEvaluated <- ExpressionEvaluator.eval(entrynameExpression, ParameterContext(inputValues), expressionLib)
        entryNameValidated <- mustBeString(entryNameExpressionEvaluated)
      } yield entryNameValidated
    }

    entry match {
      case IwdrListingArrayEntry.StringDirent(content, direntEntryName, _) => for {
        entryNameValidated <- evaluateEntryName(direntEntryName)
        writtenFile <- validate(Await.result(ioFunctionSet.writeFile(entryNameValidated, content), Duration.Inf))
      } yield writtenFile

      case IwdrListingArrayEntry.ExpressionDirent(content, direntEntryName, _) =>
        val entryEvaluation: ErrorOr[WomValue] = ExpressionEvaluator.eval(content, ParameterContext(inputValues), expressionLib)
        entryEvaluation flatMap {
          case f: WomFile =>
            val errorOrEntryName: ErrorOr[String] = direntEntryName match {
              case Some(en) => evaluateEntryName(en)
              case None => f.value.split('/').last.validNel
            }
            errorOrEntryName flatMap { entryName =>
              validate(Await.result(ioFunctionSet.copyFile(f.value, entryName), Duration.Inf))
            }
          case other => for {
            coerced <- WomStringType.coerceRawValue(other).toErrorOr
            contentString = coerced.asInstanceOf[WomString].value
            // We force the entryname to be specified, and then evaluate it:
            entryNameUnoptioned <- direntEntryName.toErrorOr("Invalid dirent: Entry was a string but no file name was supplied")
            entryname <- evaluateEntryName(entryNameUnoptioned)
            writtenFile <- validate(Await.result(ioFunctionSet.writeFile(entryname, contentString), Duration.Inf))
          } yield writtenFile
        }
      case IwdrListingArrayEntry.Expression(expression) =>
        // A single expression which must evaluate to an array of Files
        val expressionEvaluation = ExpressionEvaluator.eval(expression, ParameterContext(inputValues), expressionLib)

        expressionEvaluation flatMap {
          case array: WomArray if WomArrayType(WomSingleFileType).coercionDefined(array) => WomArrayType(WomSingleFileType).coerceRawValue(array).toErrorOr
//            val newFileArray: ErrorOr[List[WomFile]] =
//              array.value.map(_.valueString).toList.traverse[ErrorOr, WomSingleFile] { source: String =>
//                val dest = source.split('/').last
//                validate(Await.result(ioFunctionSet.copyFile(source, dest), Duration.Inf))
//              }
//            newFileArray map { nfa => WomArray(WomArrayType(WomSingleFileType), nfa) }
          case file: WomSingleFile =>
            file.validNel

          case other => s"InitialWorkDirRequirement listing expression must be Array[File] but got ${other.womType.toDisplayString}".invalidNel
        }

      case _ => ??? // TODO CWL and the rest....
    }
  }


  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
    "Programmer error: Shouldn't use InitialWorkDirRequirement listing to find output files. You silly goose.".invalidNel

  /**
    * We already get all of the task inputs when evaluating, and we don't need to highlight anything else
    */
  override def inputs: Set[String] = Set.empty
}
