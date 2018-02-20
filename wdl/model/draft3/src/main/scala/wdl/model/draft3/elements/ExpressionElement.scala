package wdl.model.draft3.elements

import wom.values.WomPrimitive

sealed trait ExpressionElement

object ExpressionElement {
  final case class PrimitiveLiteralExpressionElement(value: WomPrimitive) extends ExpressionElement

  final case class ObjectLiteral(elements: Map[String, ExpressionElement]) extends ExpressionElement
  final case class ArrayLiteral(elements: Array[ExpressionElement]) extends ExpressionElement
  final case class MapLiteral(elements: Map[String, ExpressionElement]) extends ExpressionElement
  final case class PairLiteral(left: ExpressionElement, right: ExpressionElement) extends ExpressionElement

  /**
    * Represents a unary operation (i.e. a operator symbol followed by a single argument expression)
    */
  sealed trait UnaryOperation extends ExpressionElement {
    /**
      * The expression which follows the unary operator. The argument to the operation.
      */
    def argument: ExpressionElement
  }

  final case class LogicalNot(override val argument: ExpressionElement) extends UnaryOperation
  final case class UnaryNegation(override val argument: ExpressionElement) extends UnaryOperation
  final case class UnaryBooleanNot(override val argument: ExpressionElement) extends UnaryOperation

  /**
    * A two-argument expression. Almost certainly comes from an infix operation in WDL (eg the '+' in  '7 + read_int(x)')
    */
  sealed trait BinaryOperation extends ExpressionElement {
    /**
      * The left-hand-side of the operation ('7' in the example above).
      */
    def left: ExpressionElement

    /**
      * The right-hand-side of the operation ('read_int(x)' in the example above).
      */
    def right: ExpressionElement
  }

  final case class LogicalOr(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class LogicalAnd(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class Equals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class NotEquals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class LessThan(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class LessThanOrEquals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class GreaterThan(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class GreaterThanOrEquals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class Add(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class Subtract(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class Multiply(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class Divide(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  final case class Remainder(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation

  sealed trait FunctionCall extends ExpressionElement
  final case object StdoutCall extends FunctionCall
  final case object StderrCall extends FunctionCall
  // TODO: and other engine functions

  /**
    * A single identifier lookup expression, eg Int x = y
    */
  final case class IdentifierLookup(identifier: String)

  /**
    * Represents a member access.
    *
    * But, why the split into firstIdentifier, secondIdentifierOrFirstMemberAccess and memberAccessTail? Take a look at some examples:
    *
    * ---
    *
    * Example 1: task output lookup expression:
    * Pair[Pair[String, Int], Int] pair_of_pairs = my_task.pair_of_pairs
    *
    * The first identifier is 'my_task', the second identifier is 'pair_of_pairs' and the tail is [ ]
    *
    * ---
    *
    * Example 2:
    * Int x = pair_of_pairs.left.right
    *
    * The first identifier is 'pair_of_pairs', the second is 'left' and the tail is [ 'right' ]
    *
    * ---
    *
    * Example 3:
    * Int x = mytask.pair_of_pair_output.left.right
    *
    * The first identifier is mytask, the second is pair_of_pair_output, and the tail is ['left', 'right']
    *
    * ---
    *
    * Observations:
    *  - We always get at least two identifiers.
    *  - The firstIdentifier is *always* part of the identifier we'll need to look up in the value store.
    *  - The tail will *always* be a chain of member accesses on a WomValue.
    *  - But, the second element might be part of the identifier to look up (eg my_task.pair_of_pairs) OR it might
    *      be part of a member access chain (eg pair_of_pairs.left.right). We won't know until we do the linking.
    */
  final case class IdentifierMemberAccess(firstIdentifier: String, secondIdentifierOrFirstMemberAccess: String, memberAccessTail: Vector[String]) extends ExpressionElement

  /**
    * A member access which is based on an expression rather than an identifier.
    *
    * eg:
    * (1, 2).left
    */
  final case class ExpressionMemberAccess(expression: ExpressionElement, memberAccessTail: Vector[String]) extends ExpressionElement
}
