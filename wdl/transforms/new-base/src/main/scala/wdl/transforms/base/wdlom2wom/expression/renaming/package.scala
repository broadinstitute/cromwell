package wdl.transforms.base.wdlom2wom.expression

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.transforms.base.wdlom2wom.expression.renaming.IdentifierLookupRenamer.ops._

import wdl.transforms.base.wdlom2wom.expression.renaming.BinaryOperatorEvaluators._
import wdl.transforms.base.wdlom2wom.expression.renaming.EngineFunctionEvaluators._
import wdl.transforms.base.wdlom2wom.expression.renaming.LiteralEvaluators._
import wdl.transforms.base.wdlom2wom.expression.renaming.LookupEvaluators._
import wdl.transforms.base.wdlom2wom.expression.renaming.TernaryIfEvaluator._
import wdl.transforms.base.wdlom2wom.expression.renaming.UnaryOperatorEvaluators._

package object renaming {

  implicit val expressionEvaluator: IdentifierLookupRenamer[ExpressionElement] =
    new IdentifierLookupRenamer[ExpressionElement] {
      override def renameIdentifiers(a: ExpressionElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement]
      ): ExpressionElement =
        a match {
          // Literals:
          case a: PrimitiveLiteralExpressionElement => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: StringLiteral => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: StringExpression => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ObjectLiteral => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: MapLiteral => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ArrayLiteral => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: PairLiteral => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          // Lookups and member accesses:
          case a: IdentifierLookup => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ExpressionMemberAccess => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: IdentifierMemberAccess => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: IndexAccess => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          // Unary operators:
          case a: UnaryNegation => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: UnaryPlus => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: LogicalNot => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          // Binary operators (at some point we might want to split these into separate cases):
          case a: LogicalOr => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: LogicalAnd => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Equals => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: NotEquals => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: LessThan => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: LessThanOrEquals => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: GreaterThan => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: GreaterThanOrEquals => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Add => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Subtract => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Multiply => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Divide => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Remainder => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          case a: TernaryIf => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          // Engine functions:
          case a: StdoutElement.type => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: StderrElement.type => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          case a: ReadLines => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadTsv => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadMap => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadObject => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadObjects => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadJson => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadInt => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadString => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadFloat => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: ReadBoolean => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: WriteLines => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: WriteTsv => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: WriteMap => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: WriteObject => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: WriteObjects => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: WriteJson => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Range => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Transpose => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Length => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Flatten => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: SelectFirst => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: SelectAll => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Defined => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Floor => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Ceil => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Round => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Glob => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          case a: Size => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Basename => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          case a: Zip => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Cross => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
          case a: Prefix => a.renameIdentifiers(renamingMap)(expressionElementRenamer)

          case a: Sub => a.renameIdentifiers(renamingMap)(expressionElementRenamer)
        }
    }
}
