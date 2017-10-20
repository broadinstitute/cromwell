package wom.types

import wom.values.WomValue

case object WomAnyType extends WomType {
  val toDisplayString: String = s"Any"

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
      /* This does throw an exception if it couldn't coerce (.get is intentional) */
      WomType.womTypeCoercionOrder.map(_.coerceRawValue(any)).find(_.isSuccess).getOrElse(
        throw new UnsupportedOperationException(s"Could not coerce $any into a WOM type")
      ).get
  }

  override final def isCoerceableFrom(otherType: WomType): Boolean = true
}
