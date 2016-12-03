package cromwell.core.simpleton

import wdl4s.values._

case class WdlValueSimpleton(simpletonKey: String, simpletonValue: WdlPrimitive)

/**
  * Encodes potentially complex `WdlValue`s as `Iterable[WdlValueSimpleton]`s using a protocol similar to that of
  * `MetadataValue`s: array elements are referenced by a square-brace enclosed index, map values are referenced by a
  * colon-prefixed key.  The `WdlValueSimpleton` protocol differs slightly from the `MetadataValue` protocol because map
  * keys are not specified within Cromwell: WDL authors can choose any map keys allowed by WDL, possibly including square
  * brace or colon characters that this protocol uses as metacharacters.  So the `WdlValueSimpleton` protocol escapes
  * any `[`, `]`, or `:` metacharacters in map keys with a leading backslash, which must then be unescaped whenever
  * `WdlValueSimpleton`s are transformed back to `WdlValue`s.
  */
object WdlValueSimpleton {

  implicit class KeyMetacharacterEscaper(val key: String) extends AnyVal {
    // The escapes are necessary on the first arguments to `replaceAll` since they're treated like regular expressions
    // and square braces are character class delimiters.  Backslashes must be escaped in both parameters.
    // Ignore the red in some of the "raw" strings, IntelliJ and GitHub don't seem to understand them.
    def escapeMeta = key.replaceAll(raw"\[", raw"\\[").replaceAll(raw"\]", raw"\\]").replaceAll(":", raw"\\:")
    def unescapeMeta = key.replaceAll(raw"\\\[", "[").replaceAll(raw"\\\]", "]").replaceAll(raw"\\:", ":")
  }

  implicit class WdlValueSimplifier(wdlValue: WdlValue) {
    def simplify(name: String): Iterable[WdlValueSimpleton] = wdlValue match {
      case prim: WdlPrimitive => List(WdlValueSimpleton(name, prim))
      case opt: WdlOptionalValue => opt.value.map(_.simplify(name)).getOrElse(Seq.empty)
      case WdlArray(_, arrayValue) => arrayValue.zipWithIndex flatMap { case (arrayItem, index) => arrayItem.simplify(s"$name[$index]") }
      case WdlMap(_, mapValue) => mapValue flatMap { case (key, value) => value.simplify(s"$name:${key.valueString.escapeMeta}") }
      case WdlPair(left, right) => left.simplify(s"$name:left") ++ right.simplify(s"$name:right")
      case wdlObject: WdlObjectLike => wdlObject.value flatMap { case (key, value) => value.simplify(s"$name:${key.escapeMeta}") }
      case other => throw new Exception(s"Cannot simplify wdl value $other of type ${other.wdlType}")
    }
  }

  implicit class WdlValuesSimplifier(wdlValues: Map[String, WdlValue]) {
    def simplify: Iterable[WdlValueSimpleton] = wdlValues flatMap { case (name, value) => value.simplify(name) }
  }

}
