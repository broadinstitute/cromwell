package wom.types

import wom.values.WomValue

import scala.util.{Success, Try}

case object WomAnyType extends WomType {
  val stableName: String = s"Any"

  /**
   * WomAnyType does behave slightly differently than the other WomTypes.
   * For something like WomInteger, the coercion function acts as follows:
   *
   * "Give me one of x different type of objects that are compatible with
   *  WomInteger, and I'll give you a WomInteger back"
   *
   * WomAnyType's coercion semantics are more like:
   *
   * "Give me anything at all, and I've got a heuristic that will return you
   *  the WomValue that's able to accept this Any value."
   *
   * So you give this coercion() function a String "foobar", and it returns
   * WomString("foobar"), you give it a JsNumber(2), it returns WomInteger(2).
   *
   * This is used when evaluating literal expressions. For example, a user
   * might do this in source code:
   *
   * <pre>
   * Map[Int, String] first = {"foo": 2, 3.14: "bar"}
   * Map[Float, String] second = {"foo": "bar"}
   * </pre>
   *
   * When we're simply parsing the expression {"foo": 2, 3.14: "bar"}, there
   * are two levels of coercion happening. First it just tries to coerce this
   * into ANY Map[K, V] type. It uses Map[Any, Any].coerce for that. The first
   * example above would fail this coercion stage. The second would pass and
   * return a Map[String, String].
   *
  * The second coercion that happens is the coercion to the target type.
   * In the example above, second starts out coerced as Map[String, String],
   * and then it is Map[Float, String].coerce is called on that map. This step
   * should fail at this stage.
   */
  override protected def coercion = {
    case womValue: WomValue => womValue
    case any: Any =>

      def foldFun(acc: Option[WomValue], nextType: WomType): Option[WomValue] = acc.orElse(nextType.coerceRawValue(any).toOption)

      /* This does throw an exception if it couldn't coerce (.get is intentional) */
      WomType.womTypeCoercionOrder.foldLeft[Option[WomValue]](None)(foldFun).getOrElse(
            throw new UnsupportedOperationException(s"Could not coerce $any into a WOM type"))
//
//
//      WomType.womTypeCoercionOrder.map(_.coerceRawValue(any)).find(_.isSuccess).getOrElse(
//        throw new UnsupportedOperationException(s"Could not coerce $any into a WOM type")
//      ).get
  }

  override final def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = true

  override def add(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def subtract(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def multiply(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def divide(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def mod(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def equalsType(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def notEquals(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def lessThan(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def lessThanOrEqual(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def greaterThan(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def greaterThanOrEqual(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def or(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def and(rhs: WomType): Try[WomType] = Success(WomAnyType)
  override def not: Try[WomType] = Success(WomAnyType)
  override def unaryPlus: Try[WomType] = Success(WomAnyType)
  override def unaryMinus: Try[WomType] = Success(WomAnyType)
}
