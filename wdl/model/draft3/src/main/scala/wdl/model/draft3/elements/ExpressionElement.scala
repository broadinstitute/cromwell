package wdl.model.draft3.elements

import cats.data.NonEmptyList
import wom.values.WomPrimitive

trait ExpressionElement

object ExpressionElement {
  final case class PrimitiveLiteralExpressionElement(value: WomPrimitive) extends ExpressionElement

  case object NoneLiteralElement extends ExpressionElement

  final case class StringExpression(pieces: Seq[StringPiece]) extends ExpressionElement
  sealed trait StringPiece
  final case class StringLiteral(value: String) extends StringPiece with ExpressionElement

  /**
    * For use within a StringExpression. Cannot be a standalone ExpressionElement.
    */
  final case class StringPlaceholder(expr: ExpressionElement) extends StringPiece

  sealed trait StringEscapeSequence extends StringPiece {
    def unescape: String
  }
  case object NewlineEscape extends StringEscapeSequence { override val unescape: String = System.lineSeparator }
  case object TabEscape extends StringEscapeSequence { override val unescape: String = "\t" }
  case object DoubleQuoteEscape extends StringEscapeSequence { override val unescape: String = "\"" }
  case object SingleQuoteEscape extends StringEscapeSequence { override val unescape: String = "'" }
  case object BackslashEscape extends StringEscapeSequence { override val unescape: String = "\\" }
  final case class UnicodeCharacterEscape(codePoint: Int) extends StringEscapeSequence {
    override val unescape: String = codePoint.toChar.toString
  }


  final case class KvPair(key: String, value: ExpressionElement)
  final case class ObjectLiteral(elements: Map[String, ExpressionElement]) extends ExpressionElement
  final case class ArrayLiteral(elements: Seq[ExpressionElement]) extends ExpressionElement
  final case class MapLiteral(elements: Map[ExpressionElement, ExpressionElement]) extends ExpressionElement
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
  final case class UnaryPlus(override val argument: ExpressionElement) extends UnaryOperation

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

  final case class TernaryIf(condition: ExpressionElement, ifTrue: ExpressionElement, ifFalse: ExpressionElement) extends ExpressionElement

  sealed trait FunctionCallElement extends ExpressionElement
  // 0-param functions
  case object StdoutElement extends FunctionCallElement
  case object StderrElement extends FunctionCallElement

  // 1-param functions
  sealed trait OneParamFunctionCallElement extends FunctionCallElement { def param: ExpressionElement }
  final case class Keys(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class AsMap(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class AsPairs(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class CollectByKey(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadLines(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadTsv(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadMap(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadObject(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadObjects(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadJson(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadInt(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadString(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadFloat(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class ReadBoolean(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class WriteLines(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class WriteTsv(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class WriteMap(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class WriteObject(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class WriteObjects(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class WriteJson(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Range(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Transpose(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Length(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Flatten(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class SelectFirst(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class SelectAll(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Defined(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Floor(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Ceil(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Round(param: ExpressionElement) extends OneParamFunctionCallElement
  final case class Glob(param: ExpressionElement) extends OneParamFunctionCallElement

  // 1- or 2-param functions:
  sealed trait OneOrTwoParamFunctionCallElement extends FunctionCallElement {
    def firstParam: ExpressionElement
    def secondParam: Option[ExpressionElement]
  }
  final case class Size(file: ExpressionElement, unit: Option[ExpressionElement]) extends OneOrTwoParamFunctionCallElement {
    override def firstParam: ExpressionElement = file
    override def secondParam: Option[ExpressionElement] = unit
  }
  final case class Basename(param: ExpressionElement, suffixToRemove: Option[ExpressionElement]) extends OneOrTwoParamFunctionCallElement {
    override def firstParam: ExpressionElement = param
    override def secondParam: Option[ExpressionElement] = suffixToRemove
  }

  // 2-param functions:
  sealed trait TwoParamFunctionCallElement extends FunctionCallElement {
    def arg1: ExpressionElement
    def arg2: ExpressionElement
  }
  final case class Zip(arg1: ExpressionElement, arg2: ExpressionElement) extends TwoParamFunctionCallElement
  final case class Cross(arg1: ExpressionElement, arg2: ExpressionElement) extends TwoParamFunctionCallElement
  final case class Prefix(prefix: ExpressionElement, array: ExpressionElement) extends TwoParamFunctionCallElement {
    override def arg1: ExpressionElement = prefix
    override def arg2: ExpressionElement = array
  }

  // 3-param functions:
  sealed trait ThreeParamFunctionCallElement extends FunctionCallElement {
    def arg1: ExpressionElement
    def arg2: ExpressionElement
    def arg3: ExpressionElement
  }
  final case class Sub(input: ExpressionElement, pattern: ExpressionElement, replace: ExpressionElement) extends ThreeParamFunctionCallElement {
    override def arg1: ExpressionElement = input
    override def arg2: ExpressionElement = pattern
    override def arg3: ExpressionElement = replace
  }

  /**
    * A single identifier lookup expression, eg Int x = y
    */
  final case class IdentifierLookup(identifier: String) extends ExpressionElement

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
  final case class IdentifierMemberAccess(first: String, second: String, memberAccessTail: Seq[String]) extends ExpressionElement

  /**
    * A member access which is based on an expression rather than an identifier.
    *
    * eg:
    * (1, 2).left
    */
  final case class ExpressionMemberAccess(expression: ExpressionElement, memberAccessTail: NonEmptyList[String]) extends ExpressionElement

  /**
    *
    * @param expressionElement The element being accessed
    * @param index The index expression
    */
  final case class IndexAccess(expressionElement: ExpressionElement, index: ExpressionElement) extends ExpressionElement

  /**
    * Synthetic element - never constructed from parsing an actual WDL.
    *
    * Functions as a passthrough to represent programmatically-generated WDL expressions.
    */
  final case class ExpressionLiteralElement(expression: String) extends ExpressionElement
}
