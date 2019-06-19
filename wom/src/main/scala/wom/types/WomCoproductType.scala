package wom.types
import cats.data.NonEmptyList
import wom.WomExpressionException
import wom.values.{WomCoproductValue, WomValue}

import scala.util.{Failure, Success, Try}
import mouse.all._

/**
  * Handles the possibility that a value could be one of the specified types.  At present this is only supported by CWL.
  */
case class WomCoproductType(types: NonEmptyList[WomType]) extends WomType {

  /**
    * Method to be overridden by implementation classes defining a partial function
    * for the conversion of raw input values to specific implementation class value types.
    * i.e.  `WomBooleanType` should define a partial function that knows how to
    * construct `WomBoolean`s for inputs of supported types and contents.  Values for which
    * the partial function is not defined are assumed to not be convertible to the target type.
    */
  override def coercion: PartialFunction[Any, WomValue] = {
    case wct@WomCoproductValue(tpe, _) if (tpe.equalsType(this).isSuccess) => wct

    //If we can find this type exactly in our coproduct, use that type for the coercion
    case womValue: WomValue if (types.toList.contains(womValue.womType)) =>
      val v: Try[WomValue] = womValue.womType.coerceRawValue(womValue)
      WomCoproductValue(this, v.get)

    //If we don't have any information, try to coerce this value one by one, stopping at the first successful try
    case any =>
      val triedToCoerce: Try[WomValue] = types.map(_.coerceRawValue(any)).toList.reduce(_ orElse _)

      triedToCoerce.getOrElse(throw new WomTypeException(s"unable to coerce $any to a member of the set of types ${types.toList.mkString(", ")}")) |> (WomCoproductValue(this, _))
  }

  def typeExists(tpe: WomType): Try[WomBooleanType.type] =
    types.exists(_.equals(tpe)) match {
      case true => Success(WomBooleanType)
      case _ => Failure(new WomExpressionException(s"Type equality assertion failed because $tpe was not found in the coproduct of ${stableName}"))
    }

  override def stableName: String =
    types.map(_.stableName).toList.mkString("Coproduct[",", ", "]")

  override def equalsType(rhs: WomType): Try[WomType] =
    rhs match {
      case WomCoproductType(tpes) if types.equals(tpes) => Success(WomBooleanType)
      case _ if types.exists(_.equals(rhs)) =>  Success(WomBooleanType)
      case _ => Failure(new WomExpressionException(s"Type equality could not be asserted because $rhs was found in the coproduct of $stableName"))
    }


}
