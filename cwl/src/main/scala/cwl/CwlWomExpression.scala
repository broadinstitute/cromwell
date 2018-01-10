package cwl

import cats.data.NonEmptyList
import cats.syntax.either._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.Validation._
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import mouse.all._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.util.Try

trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WomType

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = cwlExpressionType.validNel
}

case class JobPreparationExpression(expression: Expression,
                                    override val inputs: Set[String]) extends CwlWomExpression {
  val cwlExpressionType = WomAnyType

  override def sourceString = expression match {
    case Expression.ECMAScriptExpression(s) => s.value
    case Expression.ECMAScriptFunction(s) => s.value
  }

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    val pc = ParameterContext(inputValues)
    expression.
      fold(EvaluateExpression).
      apply(pc).
      cata(Right(_), Left(_)). // this is because toEither is not a thing in scala 2.11.
      leftMap(e => NonEmptyList.one(e.getMessage)).
      toValidated
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = Set.empty[WomFile].validNel
}

final case class InitialWorkDirFileGeneratorExpression(entry: IwdrListingArrayEntry) extends CwlWomExpression {
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
        entryNameExpressionEvaluated <- ExpressionEvaluator.evalCwlExpression(entrynameExpression, ParameterContext(inputValues)).toErrorOr
        entryNameValidated <- mustBeString(entryNameExpressionEvaluated)
      } yield entryNameValidated
    }

    entry match {
      case IwdrListingArrayEntry.StringDirent(content, direntEntryName, _) => for {
        entryNameValidated <- evaluateEntryName(direntEntryName)
        writtenFile <- Try(Await.result(ioFunctionSet.writeFile(entryNameValidated, content), Duration.Inf)).toErrorOr
      } yield writtenFile

      case IwdrListingArrayEntry.ExpressionDirent(contentExpression, direntEntryName, _) =>
        val entryEvaluation: ErrorOr[WomValue] = ExpressionEvaluator.evalCwlExpression(contentExpression, ParameterContext().addInputs(inputValues)).toErrorOr
        entryEvaluation flatMap {
          case f: WomFile =>
            val entryName: ErrorOr[String] = direntEntryName match {
              case Some(en) => evaluateEntryName(en)
              case None => f.value.split('/').last.validNel
            }
            entryName flatMap { en => Try(Await.result(ioFunctionSet.copyFile(f.value, en), Duration.Inf)).toErrorOr }
          case other => for {
            coerced <- WomStringType.coerceRawValue(other).toErrorOr
            contentString = coerced.asInstanceOf[WomString].value
            // We force the entryname to be specified, and then evaluate it:
            entryNameUnoptioned <- direntEntryName.toErrorOr("Invalid dirent: Entry was a string but no file name was supplied")
            entryname <- evaluateEntryName(entryNameUnoptioned)
            writtenFile <- Try(Await.result(ioFunctionSet.writeFile(entryname, contentString), Duration.Inf)).toErrorOr
          }  yield writtenFile
        }
      case IwdrListingArrayEntry.Expression(expression) =>
        // A single expression which must evaluate to an array of Files
        val expressionEvaluation: ErrorOr[WomValue] = ExpressionEvaluator.evalCwlExpression(expression, ParameterContext().addInputs(inputValues)).toErrorOr

        expressionEvaluation flatMap {
          case array: WomArray if WomArrayType(WomSingleFileType).coercionDefined(array) =>
            val newFileArray: ErrorOr[List[WomFile]] = array.value.map(_.valueString).toList.traverse { source: String =>
              val dest = source.split('/').last
              Try(Await.result(ioFunctionSet.copyFile(source, dest), Duration.Inf)).toErrorOr
            }
            newFileArray map { nfa => WomArray(WomArrayType(WomSingleFileType), nfa) }

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
