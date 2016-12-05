package wdl4s.types

import wdl4s.AstTools.EnhancedAstNode
import wdl4s.values._
import wdl4s.{WdlExpressionException, WdlSource, WdlSyntaxErrorFormatter}
import wdl4s.parser.WdlParser

import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}

class WdlTypeException(message: String) extends RuntimeException(message)

trait WdlType {

  /**
   * Method to be overridden by implementation classes defining a partial function
   * for the conversion of raw input values to specific implementation class value types.
   * i.e.  `WdlBooleanType` should define a partial function that knows how to
   * construct `WdlBoolean`s for inputs of supported types and contents.  Values for which
   * the partial function is not defined are assumed to not be convertible to the target type.
   */
  protected def coercion: PartialFunction[Any, WdlValue]
  def coercionDefined(any: Any) = coercion.isDefinedAt(any)

  /**
   * Public interface for a `Try`-wrapped conversion of an input of type `Any` to
   * a `WdlValue`.
   */
  def coerceRawValue(any: Any): Try[WdlValue] = {
    any match {
      case v: WdlValue if v.wdlType == this => Success(v)
      case a if !coercion.isDefinedAt(any) => Failure(new IllegalArgumentException(s"No coercion defined from '$any' of type '${any.getClass}' to ${getClass.getSimpleName}."))
      case _ => Try(coercion(any))
    }
  }

  def isCoerceableFrom(otherType: WdlType): Boolean = false

  def toWdlString: String

  /**
   * Converts WDL source into a WdlValue of this type, if possible.
   *
   * @param wdlSource source code representing the WdlValue
   * @return The WdlValue
   */
  //TODO: return a Try ?
  def fromWdlString(wdlSource: WdlSource): WdlValue = {
    val tokens = WdlType.parser.lex(wdlSource, "string")
    val terminalMap = tokens.asScala.toVector.map {(_, wdlSource)}.toMap
    val wdlSyntaxErrorFormatter = new WdlSyntaxErrorFormatter(terminalMap)

    /* Parsing as an expression is not sufficient... only a subset of these
     * ASTs are valid as WdlValues and this distinction is done in the
     * .wdlValue() method.
     */
    val ast = WdlType.parser.parse_e(tokens, wdlSyntaxErrorFormatter).toAst

    ast.wdlValue(this, wdlSyntaxErrorFormatter)
  }

  def invalid(operation: String) = Failure(new WdlExpressionException(s"Cannot perform operation: $operation"))
  def add(rhs: WdlType): Try[WdlType] = invalid(s"$this + $rhs")
  def subtract(rhs: WdlType): Try[WdlType] = invalid(s"$this - $rhs")
  def multiply(rhs: WdlType): Try[WdlType] = invalid(s"$this * $rhs")
  def divide(rhs: WdlType): Try[WdlType] = invalid(s"$this / $rhs")
  def mod(rhs: WdlType): Try[WdlType] = invalid(s"$this % $rhs")
  def equals(rhs: WdlType): Try[WdlType] = invalid(s"$this == $rhs")
  def notEquals(rhs: WdlType): Try[WdlType] = equals(rhs) map {x => WdlBooleanType}
  def lessThan(rhs: WdlType): Try[WdlType] = invalid(s"$this < $rhs")
  def lessThanOrEqual(rhs: WdlType): Try[WdlType] = (lessThan(rhs), equals(rhs)) match {
    case (Success(b:WdlType), _) if b == WdlBooleanType => Success(WdlBooleanType)
    case (_, Success(b:WdlType)) if b == WdlBooleanType => Success(WdlBooleanType)
    case (_, _) => invalid(s"$this <= $rhs")
  }
  def greaterThan(rhs: WdlType): Try[WdlType] = invalid(s"$this > $rhs")
  def greaterThanOrEqual(rhs: WdlType): Try[WdlType] = (greaterThan(rhs), equals(rhs)) match {
    case (Success(b:WdlType), _) if b == WdlBooleanType => Success(WdlBooleanType)
    case (_, Success(b:WdlType)) if b == WdlBooleanType => Success(WdlBooleanType)
    case (_, _) => invalid(s"$this >= $rhs")
  }
  def or(rhs: WdlType): Try[WdlType] = invalid(s"$this || $rhs")
  def and(rhs: WdlType): Try[WdlType] = invalid(s"$this && $rhs")
  def not: Try[WdlType] = invalid(s"!$this")
  def unaryPlus: Try[WdlType] = invalid(s"+$this")
  def unaryMinus: Try[WdlType] = invalid(s"-$this")
}

object WdlType {
  val parser = new WdlParser()

  /* This is in the order of coercion from non-wdl types */
  val wdlTypeCoercionOrder: Seq[WdlType] = Seq(
    WdlStringType, WdlIntegerType, WdlFloatType, WdlMapType(WdlAnyType, WdlAnyType),
    WdlArrayType(WdlAnyType), WdlBooleanType, WdlObjectType
  )

  def homogeneousTypeFromValues(values: Iterable[WdlValue]): WdlType =
    homogeneousTypeFromTypes(values.map(_.wdlType))

  def homogeneousTypeFromTypes(types: Iterable[WdlType]): WdlType = {
    types.toSet match {
      case s if s.isEmpty => WdlAnyType
      case s if s.size == 1 => s.head
      case _ => lowestCommonSubtype(types)
    }
  }

  def lowestCommonSubtype(types: Iterable[WdlType]): WdlType = {
    types.collectFirst {
      case t1 if types.forall(t2 => t1.isCoerceableFrom(t2)) => t1
    } getOrElse WdlAnyType
  }

  def fromWdlString(wdlString: String): WdlType = {
    wdlString match {
      case "Expression" => WdlExpressionType
      case _ =>
        val tokens = parser.lex(wdlString, "string")
        val terminalMap = tokens.asScala.toVector.map {(_, wdlString)}.toMap
        val wdlSyntaxErrorFormatter = new WdlSyntaxErrorFormatter(terminalMap)

        /* parse_type_e() is the parse function for the $type_e nonterminal in grammar.hgr */
        val ast = parser.parse_type_e(tokens, wdlSyntaxErrorFormatter).toAst

        ast.wdlType(wdlSyntaxErrorFormatter)
    }
  }
}
