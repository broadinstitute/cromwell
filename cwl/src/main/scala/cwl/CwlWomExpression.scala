package cwl

import cats.syntax.option._
import cats.data.NonEmptyList
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import common.validation.ErrorOr.ShortCircuitingFlatMap
import cats.syntax.validated._
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.util.Try
import cats.syntax.either._
import mouse.all._

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

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) =
    ParameterContext().
      addInputs(inputValues).
      flatMap(pc =>
        expression.
          fold(EvaluateExpression).
            apply(pc).
            cata(Right(_),Left(_)). // this is because toEither is not a thing in scala 2.11.
            leftMap(e => NonEmptyList.one(e.getMessage))
      ).toValidated

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = Set.empty[WomFile].validNel
}

final case class InitialWorkDirFileGeneratorExpression(entry: IwdrListingArrayEntry) extends CwlWomExpression {
  override def cwlExpressionType: WomType = WomSingleFileType
  override def sourceString: String = entry.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomFile] = {
    def mustBeString(womValue: WomValue): ErrorOr[String] = womValue match {
      case WomString(s) => s.validNel
      case other => WomStringType.coerceRawValue(other).map(_.asInstanceOf[WomString].value).toErrorOr
    }

    def evaluateEntryName(stringOrExpression: StringOrExpression): ErrorOr[String] = stringOrExpression match {
      case StringOrExpression.String(s) => s.validNel
      case StringOrExpression.ECMAScriptExpression(entrynameExpression) => for {
        entryNameExpressionEvaluated <- ExpressionEvaluator.evalExpression(entrynameExpression, ParameterContext().withInputs(inputValues, ioFunctionSet)).toErrorOr
        entryNameValidated <- mustBeString(entryNameExpressionEvaluated)
      } yield entryNameValidated
    }

    entry match {
      case IwdrListingArrayEntry.StringDirent(content, direntEntryName, _) => for {
        entryNameValidated <- evaluateEntryName(direntEntryName)
        writtenFile <- Try(Await.result(ioFunctionSet.writeFile(entryNameValidated, content), Duration.Inf)).toErrorOr
      } yield writtenFile

      case IwdrListingArrayEntry.ExpressionDirent(Expression.ECMAScriptExpression(contentExpression), direntEntryName, _) =>
        val entryEvaluation: ErrorOr[WomValue] = ExpressionEvaluator.evalExpression(contentExpression, ParameterContext().withInputs(inputValues, ioFunctionSet)).toErrorOr
        entryEvaluation flatMap {
          // TODO CWL: Once files have "local paths", we will be able to specify a new local name based on direntEntryName if necessary.
          case f: WomFile => f.validNel
          case other => for {
            coerced <- WomStringType.coerceRawValue(other).toErrorOr
            contentString = coerced.asInstanceOf[WomString].value
            // We force the entryname to be specified, and then evaluate it:
            entryNameUnoptioned <- direntEntryName.toErrorOr("Invalid dirent: Entry was a string but no file name was supplied")
            entryname <- evaluateEntryName(entryNameUnoptioned)
            writtenFile <- Try(Await.result(ioFunctionSet.writeFile(entryname, contentString), Duration.Inf)).toErrorOr
          }  yield writtenFile
        }

      case _ => ??? // TODO WOM and the rest....
    }
  }


  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
    "Programmer error: Shouldn't use InitialWorkDirRequirement listing to find output files. You silly goose.".invalidNel

  /**
    * We already get all of the task inputs when evaluating, and we don't need to highlight anything else
    */
  override def inputs: Set[String] = Set.empty
}
