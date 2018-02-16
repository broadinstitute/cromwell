package wdl.model.draft3.elements

import wom.values.WomPrimitive

sealed trait ExpressionElement

object ExpressionElement {
  case class PrimitiveLiteralExpressionElement(value: WomPrimitive) extends ExpressionElement

  case class ObjectLiteral(elements: Map[String, ExpressionElement])
  case class ArrayLiteral(elements: Array[ExpressionElement])
  case class MapLiteral(elements: Array[ExpressionElement])
  case class TupleLiteral(elements: Array[ExpressionElement])

  /**
    * Represents a unary operation (i.e. a operator symbol followed by a single argument expression)
    */
  sealed trait UnaryOperation {
    /**
      * The expression which follows the unary operator. The argument to the operation.
      */
    def argument: ExpressionElement
  }

  case class LogicalNot(override val argument: ExpressionElement) extends UnaryOperation
  case class UnaryNegation(override val argument: ExpressionElement) extends UnaryOperation
  case class UnaryBooleanNot(override val argument: ExpressionElement) extends UnaryOperation

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

  case class LogicalOr(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class LogicalAnd(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class Equals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class NotEquals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class LessThan(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class LessThanOrEquals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class GreaterThan(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class GreaterThanOrEquals(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class Add(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class Subtract(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class Multiply(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class Divide(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation
  case class Remainder(override val left: ExpressionElement, override val right: ExpressionElement) extends BinaryOperation

  sealed trait FunctionCall extends ExpressionElement
  case object StdoutCall extends FunctionCall
  case object StderrCall extends FunctionCall
  // TODO: and the rest...



}