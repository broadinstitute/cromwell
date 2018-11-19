package wdl.transforms.biscayne.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.transforms.base.linking.expression.files.BinaryOperatorEvaluators._
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators._
import wdl.transforms.base.linking.expression.files.LiteralEvaluators._
import wdl.transforms.base.linking.expression.files.LookupEvaluators._
import wdl.transforms.base.linking.expression.files.TernaryIfEvaluator.ternaryIfEvaluator
import wdl.transforms.base.linking.expression.files.UnaryOperatorEvaluators._
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}
import wdl.transforms.biscayne.linking.expression.files.BiscayneFileEvaluators._

package object files {

  implicit val expressionFileEvaluator: FileEvaluator[ExpressionElement] = new FileEvaluator[ExpressionElement] {

    override def predictFilesNeededToEvaluate(a: ExpressionElement, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {

      a match {
        // Literals:
        case a: PrimitiveLiteralExpressionElement => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: StringLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ObjectLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: MapLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ArrayLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: PairLiteral => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        // Lookups and member accesses:
        case a: IdentifierLookup => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ExpressionMemberAccess => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: IdentifierMemberAccess => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: IndexAccess => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        // Unary operators:
        case a: UnaryNegation => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: UnaryPlus => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: LogicalNot => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        // Binary operators (at some point we might want to split these into separate cases):
        case a: LogicalOr => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: LogicalAnd => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Equals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: NotEquals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: LessThan => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: LessThanOrEquals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: GreaterThan => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: GreaterThanOrEquals => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Add => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Subtract => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Multiply => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Divide => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Remainder => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        case a: TernaryIf => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        // Engine functions:
        case StdoutElement => StdoutElement.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case StderrElement => StderrElement.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        case a: ReadLines => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadTsv => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadMap => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadObject => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadObjects => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadJson => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadInt => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadString => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadFloat => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: ReadBoolean => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: WriteLines => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: WriteTsv => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: WriteMap => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: WriteObject => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: WriteObjects => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: WriteJson => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Range => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Transpose => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Length => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Flatten => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Prefix => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: SelectFirst => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: SelectAll => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Defined => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Floor => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Ceil => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Round => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Glob => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        case a: Size => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Basename => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        case a: Zip => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: Cross => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        case a: Sub => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        case a: Keys => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: AsMap => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: AsPairs => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)
        case a: CollectByKey => a.predictFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)(fileEvaluator, valueEvaluator)

        case other => s"No implementation of FileEvaluator[${other.getClass.getSimpleName}]".invalidNel
      }
    }
  }
}
