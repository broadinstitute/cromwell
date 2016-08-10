package cromwell.core.simpleton

import wdl4s.values.{WdlArray, WdlMap, WdlObjectLike, WdlPrimitive, WdlValue}

case class WdlValueSimpleton(simpletonKey: String, simpletonValue: WdlPrimitive)

object WdlValueSimpleton {

  implicit class WdlValueSimplifier(wdlValue: WdlValue) {
    def simplify(name: String): Iterable[WdlValueSimpleton] = wdlValue match {
      case prim: WdlPrimitive => List(WdlValueSimpleton(name, prim))
      case WdlArray(_, arrayValue) => arrayValue.zipWithIndex flatMap { case (arrayItem, index) => arrayItem.simplify(s"$name[$index]") }
      case WdlMap(_, mapValue) => mapValue flatMap { case (key, value) => value.simplify(s"$name:${key.valueString}") }
      case wdlObject: WdlObjectLike => wdlObject.value flatMap { case (key, value) => value.simplify(s"$name:$key") }
    }
  }

  implicit class WdlValuesSimplifier(wdlValues: Map[String, WdlValue]) {
    def simplify: Iterable[WdlValueSimpleton] = wdlValues flatMap { case (name, value) => value.simplify(name) }
  }
}
