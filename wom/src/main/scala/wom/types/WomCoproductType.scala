package wom.types
import common.validation.ErrorOr.ErrorOr
import wom.WomExpressionException
import wom.values.{WomOptionalValue, WomValue}

import scala.util.{Failure, Success, Try}

case class WomCoproductType(types: Set[WomType]) extends WomType {
  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WomBooleanType` should define a partial function that knows how to
    * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override def coercion: PartialFunction[Any, WomValue] = {
    case WomOptionalValue(tpe, Some(value)) =>
      //todo: Option . get here but I dont know how to override PartialFunction's isDefined to make it not evaluate if type isn't present
      //todo: Try.get for same reasons
      types.find(_ == tpe).get.coerceRawValue(value).get
    case any =>
      val f: PartialFunction[Any, WomValue] = types.map(
        t =>
          PartialFunction.apply[Any, WomValue](any =>
            t.coerceRawValue(any).get)
      ).reduce(_ orElse _)

      f(any)
  }

  def typeExists(tpe: WomType): Try[WomBooleanType.type] =
    types.exists(_.equals(WomStringType)) match {
      case true => Success(WomBooleanType)
      case _ => Failure(new WomExpressionException(s"Type equality could not be asserted because $tpe was not found in the coproduct of ${toDisplayString}"))
    }

  override def toDisplayString: String =
    types.map(_.toDisplayString).mkString("Coproduct[",", ", "]")

  override def equals(rhs: WomType): Try[WomType] =
    types.exists(_.equals(rhs)) match {
      case true => Success(WomBooleanType)
      case _ => Failure(new WomExpressionException(s"Type equality could not be asserted because $rhs was found in the coproduct of $toDisplayString"))
    }
}
