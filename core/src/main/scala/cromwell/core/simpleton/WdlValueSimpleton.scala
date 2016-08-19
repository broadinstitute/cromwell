package cromwell.core.simpleton

import wdl4s.values.{WdlArray, WdlMap, WdlObjectLike, WdlPrimitive, WdlValue}

case class WdlValueSimpleton(simpletonKey: String, simpletonValue: WdlPrimitive)

object WdlValueSimpleton {

  implicit class KeyMetacharacterEscaper(val key: String) extends AnyVal {
    // The escapes are necessary on the first arguments to `replaceAll` since they're treated like regular expressions
    // and square braces are character class delimiters.  Backslashes must be escaped in both parameters.
    // Ignore the red in some of the "raw" strings, Intellij and GitHub don't seem to understand them.
    def escapeMeta = key.replaceAll(raw"\[", raw"\\[").replaceAll(raw"\]", raw"\\]").replaceAll(":", raw"\\:")
    def unescapeMeta = key.replaceAll(raw"\\\[", "[").replaceAll(raw"\\\]", "]").replaceAll(raw"\\:", ":")
  }

  implicit class WdlValueSimplifier(wdlValue: WdlValue) {
    def simplify(name: String): Iterable[WdlValueSimpleton] = wdlValue match {
      case prim: WdlPrimitive => List(WdlValueSimpleton(name, prim))
      case WdlArray(_, arrayValue) => arrayValue.zipWithIndex flatMap { case (arrayItem, index) => arrayItem.simplify(s"$name[$index]") }
      case WdlMap(_, mapValue) => mapValue flatMap { case (key, value) => value.simplify(s"$name:${key.valueString.escapeMeta}") }
      case wdlObject: WdlObjectLike => wdlObject.value flatMap { case (key, value) => value.simplify(s"$name:${key.escapeMeta}") }
    }
  }

  implicit class WdlValuesSimplifier(wdlValues: Map[String, WdlValue]) {
    def simplify: Iterable[WdlValueSimpleton] = wdlValues flatMap { case (name, value) => value.simplify(name) }
  }

}
