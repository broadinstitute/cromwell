package wom.util

import jdk.nashorn.api.scripting.AbstractJSObject

import scala.collection.JavaConverters._

/**
  * Builds a JSObject compatible wrapper around an array.
  *
  * Without this wrapper, JSON.stringify may not work, as the native java array does not implement getMember.
  */
final case class JsArray(seq: Seq[AnyRef]) extends AbstractJSObject {
  override def isArray: Boolean = true

  override def getSlot(index: Int): AnyRef = seq(index)

  override def hasSlot(slot: Int): Boolean = seq.lengthCompare(slot) < 0

  override def values(): java.util.Collection[AnyRef] = seq.asJava

  override def getMember(name: String): AnyRef = {
    name match {
      case "length" => Int.box(seq.size)
      case _ => super.getMember(name)
    }
  }
}
