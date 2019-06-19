package wom.types

import spray.json._
import wom.values.{WomPair, WomValue}

import scala.util.{Failure, Success, Try}

case class WomPairType(leftType: WomType, rightType: WomType) extends WomType {

  override def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = otherType match {
    case WomPairType(otherType1, otherType2) => leftType.isCoerceableFrom(otherType1) && rightType.isCoerceableFrom(otherType2)
    case _ => false
  }

  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WomBooleanType` should define a partial function that knows how to
    * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override protected def coercion: PartialFunction[Any, WomValue] = {
    case otherPair @ WomPair(otherValue1, otherValue2) if isCoerceableFrom(otherPair.womType) =>
      WomPair(leftType.coerceRawValue(otherValue1).get, rightType.coerceRawValue(otherValue2).get)
    case jsObject: JsObject if jsObject.fields.size == 2 => coercePair(jsObject.fields, this)
  }

  def coercePair(m: Map[String, JsValue], womPairType: WomPairType): WomPair = {

    val caseNormalizedMap = m map { case(k, v) => k.toLowerCase.capitalize -> v }

    def invalidPair(missingArgs: String*) = Seq(Failure(new IllegalArgumentException(s"Pair ${JsObject(m)} requires for ${missingArgs.mkString("/")} value(s) to be defined.")))

    val womPair: Seq[Try[WomValue]] = (caseNormalizedMap.get("Left"), caseNormalizedMap.get("Right")) match {
      case (Some(leftVal), Some(rightVal)) => Seq(womPairType.leftType.coerceRawValue(leftVal), womPairType.rightType.coerceRawValue(rightVal))
      case (Some(_), _) => invalidPair("Right")
      case (_, Some(_)) => invalidPair("Left")
      case _ => invalidPair("Right", "Left")
    }

    val failures = womPair collect { case f:Failure[_] => f }

    if (failures.isEmpty) {
      womPair match {
        case Seq(Success(left), Success(right)) => WomPair(left, right)
      }
    } else {
      throw new UnsupportedOperationException(s"Failed to coerce one or more values for creating a ${womPairType.stableName}:\n${failures.toList.mkString("\n")}")
    }
  }

  override def stableName: String = s"Pair[${leftType.stableName}, ${rightType.stableName}]"
}
