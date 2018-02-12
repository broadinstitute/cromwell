package wom.types
import wom.WomExpressionException
import wom.values.WomValue

import scala.util.{Failure, Success, Try}

case class WomCoproductType(types: Set[WomType]) extends WomType {
  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WomBooleanType` should define a partial function that knows how to
    * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override def coercion(): PartialFunction[Any, WomValue] =
      types.map(
        _.coercion()
      ).reduce(_ orElse _)

  override def toDisplayString: String =
    types.map(_.toDisplayString).mkString("Coproduct[",", ", "]")

  override def equals(rhs: WomType): Try[WomType] = {
    println("running equals")
    types.exists(_.equals(rhs)) match {
      case true => Success(WomBooleanType)
      case _ => Failure(new WomExpressionException(s"Type equality could not be asserted because $rhs was found in the coproduct of $toDisplayString"))
    }
  }
}
