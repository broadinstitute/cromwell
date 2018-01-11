package wom.util

import jdk.nashorn.api.scripting.AbstractJSObject

import scala.collection.JavaConverters._

/**
  * Builds a JSObject compatible wrapper around a map.
  *
  * Without this wrapper, JSON.stringify may not work, as the NativeJSON.stringify only likes JSObject and not Map.
  *
  * NOTE: With this wrapper, the Map keys _must_ be strings. Also non-string expressions are not implicitly coerced to
  * property identifiers by Nashorn. Ex: 'myMap[true]' will NOT work, but 'myMap["true"]' will work.
  */
final case class JsMap(map: Map[String, AnyRef]) extends AbstractJSObject {
  override def keySet(): java.util.Set[String] = map.keySet.asJava

  override def values(): java.util.Collection[AnyRef] = map.values.toList.asJava

  override def hasMember(name: String): Boolean = map.contains(name)

  override def getMember(name: String): AnyRef = map.getOrElse(name, super.getMember(name))
}
