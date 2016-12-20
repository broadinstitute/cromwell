package wdl4s.types

import wdl4s.values.WdlValue

case object WdlAnyType extends WdlType {
  val toWdlString: String = s"Any"

  /**
   * WdlAnyType does behave slightly differently than the other WdlTypes.
   * For something like WdlInteger, the coercion function acts as follows:
   *
   * "Give me one of x different type of objects that are compatible with
   *  WdlInteger, and I'll give you a WdlInteger back"
   *
   * WdlAnyType's coercion semantics are more like:
   *
   * "Give me anything at all, and I've got a heuristic that will return you
   *  the WdlValue that's able to accept this Any value."
   *
   * So you give this coercion() function a String "foobar", and it returns
   * WdlString("foobar"), you give it a JsNumber(2), it returns WdlInteger(2).
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
    case wdlValue: WdlValue => wdlValue
    case any: Any =>
      /* This does throw an exception if it couldn't coerce (.get is intentional) */
      WdlType.wdlTypeCoercionOrder.map(_.coerceRawValue(any)).find(_.isSuccess).getOrElse(
        throw new UnsupportedOperationException(s"Could not coerce $any into a WDL type")
      ).get
  }

  override final def isCoerceableFrom(otherType: WdlType): Boolean = true
}
