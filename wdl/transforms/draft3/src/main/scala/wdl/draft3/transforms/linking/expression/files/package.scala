package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.values.{WomFile, WomValue}
import wdl.draft3.transforms.linking.expression.files.UnaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.files.LiteralEvaluators._
import wdl.draft3.transforms.linking.expression.files.UnaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.files.BinaryOperatorEvaluators._
import wdl.draft3.transforms.linking.expression.files.LookupEvaluators._
import wdl.draft3.transforms.linking.expression.files.TernaryIfEvaluator.ternaryIfEvaluator
import wdl.draft3.transforms.linking.expression.files.EngineFunctionEvaluators._
import wom.types.WomType

package object files {

  implicit val expressionFileEvaluator: FileEvaluator[ExpressionElement] = new FileEvaluator[ExpressionElement] {

    override def predictFilesNeededToEvaluate(a: ExpressionElement, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = {

      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: StringLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ObjectLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: MapLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ArrayLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: PairLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ExpressionMemberAccess => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: IdentifierMemberAccess => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        // Unary operators:
        case a: UnaryNegation => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: UnaryPlus => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: LogicalNot => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: LogicalAnd => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Equals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: NotEquals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: LessThan => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: LessThanOrEquals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: GreaterThan => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: GreaterThanOrEquals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Add => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Subtract => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Multiply => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Divide => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Remainder => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        case a: TernaryIf => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        // Engine functions:
        case StdoutElement => StdoutElement.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case StderrElement => StderrElement.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        case a: ReadLines => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadTsv => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadMap => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadObject => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadObjects => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadJson => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadInt => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadString => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadFloat => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: ReadBoolean => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: WriteLines => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: WriteTsv => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: WriteMap => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: WriteObject => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: WriteObjects => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: WriteJson => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Range => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Transpose => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Length => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Flatten => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Prefix => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: SelectFirst => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: SelectAll => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Defined => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Floor => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Ceil => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Round => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        case a: Size => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Basename => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        case a: Zip => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
        case a: Cross => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        case a: Sub => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)

        case other => s"No implementation of FileEvaluator[${other.getClass.getSimpleName}]".invalidNel
      }
    }
  }
}
