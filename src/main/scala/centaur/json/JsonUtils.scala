package centaur.json

import spray.json.{JsArray, JsObject, JsValue}

object JsonUtils {
  implicit class EnhancedJsValue(val jsValue: JsValue) extends AnyVal {
    /**
      * Modified from http://stackoverflow.com/a/31592156 - changes were made both to port from using
      * Play-JSON to Spray-JSON as well as to handle some cases specific to Cromwell's metadata response
      *
      * Will flatten out a JsValue into a JsObject such that the keys represent the original hierarchy, e.g.
      * foo.bar.blah = 5 would have originally looked like:
      *
      * "foo": {
      *    "bar": {
      *        "blah": 5
      *    }
      * }
      *
      * The first two cases are specific to our metadata structures in that it's trying to do something smart with
      * the calls.
      *
      * The first case is a standard call and that will itself be flattened and added to the accumulator. In this
      * situation there's an array of length one with a call packet inside. It'll be detected, extracted and flattened
      *
      * The second handles scatters, where each scatter has its own element in the array. In this case all of the
      * scatters are flattened themselves with the scatter index being inserted into the flattened key name.
      */
    def flatten(prefix: String = ""): JsObject = jsValue.asJsObject.fields.foldLeft(JsObject.empty) {
      case (acc, (k, v: JsArray)) if v.isSingleCallArray => acc ++ JsObject(k -> v.elements.head).flatten()
      case (acc, (k, v: JsArray)) if v.hasShardIndex =>
        // The .get on the shardIndex is safe as we know all elements of the array have a shardIndex field
        acc ++ v.elements.map(_.asJsObject).fold(JsObject.empty) { (x, y) => x ++ y.flatten(s"$k.${y.shardIndex.get}") }
      case (acc, (k, v: JsObject)) =>
        if (prefix.isEmpty) acc.mergeWith(v.flatten(k))
        else acc.mergeWith(v.flatten(s"$prefix.$k"))
      case (acc, (k, v)) =>
        if (prefix.isEmpty) acc + (k -> v)
        else acc + (s"$prefix.$k" -> v)
    }
  }

  implicit class EnhancedJsObject(val jsObject: JsObject) extends AnyVal {
    // mergeWith and ++ are stolen from https://github.com/spray/spray-json/pull/135/
    def mergeWith(other: JsObject) = new JsObject(jsObject.fields ++ other.fields)
    def ++(other: JsObject) = this mergeWith other
    // + is "heavily inspired" by Play-JSON
    def +(otherField: (String, JsValue)) = JsObject(jsObject.fields + (otherField._1 -> otherField._2))

    // A couple of helper functions to assist with flattening Cromwell metadata responses
    def hasShardIndex = jsObject.fields.keySet contains "shardIndex"
    def shardIndex = jsObject.fields.get("shardIndex") map { _.toString() }
    def flattenToMap: Map [String, JsValue] = jsObject.flatten().fields map { case (k, v: JsValue) => k -> v}
  }

  /**
    * A collection of helper functions to assist with flattening Cromwell metadata responses
    */
  implicit class EnhancedJsArray(val jsArray: JsArray) extends AnyVal {
    def isObjectArray: Boolean = jsArray.elements forall { _.isInstanceOf[JsObject] }
    def nonEmpty = jsArray.elements.nonEmpty
    def size = jsArray.elements.size
    def nonEmptyObjectArray = jsArray.isObjectArray && jsArray.nonEmpty
    def isSingleCallArray = jsArray.hasShardIndex && jsArray.size == 1

    def hasShardIndex: Boolean = {
      if (jsArray.nonEmptyObjectArray) jsArray.elements.map(_.asJsObject) forall { _.hasShardIndex }
      else false
    }
  }
}
