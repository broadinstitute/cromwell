package wom.types
import wom.values.WomValue

import scala.util.Try

case class WomCoproductType(types: Set[WomType]) extends WomType {
  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WomBooleanType` should define a partial function that knows how to
    * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override def coercion(): PartialFunction[Any, WomValue] = {
    case any =>
      types.map(t =>
        t.coerceRawValue(any)
      ).reduce(_ orElse _).get
  }

  override def toDisplayString: String =
    types.map(_.toDisplayString).mkString("Coproduct[",", ", "]")
}
