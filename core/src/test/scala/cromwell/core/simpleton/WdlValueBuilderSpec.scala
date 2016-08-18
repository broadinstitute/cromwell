package cromwell.core.simpleton

import cromwell.core.simpleton.WdlValueBuilderSpec._
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlMapType, WdlStringType}
import wdl4s.values.{WdlArray, WdlInteger, WdlMap, WdlString}
import wdl4s.{TaskOutput, WdlExpression}

object WdlValueBuilderSpec {
  // WdlValueBuilder doesn't care about this expression, but something needs to be passed to the TaskOutput constructor.
  val IgnoredExpression = WdlExpression.fromString(""" "" """)
}

class WdlValueBuilderSpec extends FlatSpec with Matchers {

  "Builder" should "build" in {

    val wdlValues = Map(
      "foo" -> WdlString("none"),
      "bar" -> WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(2))),
      "baz" -> WdlArray(WdlArrayType(WdlArrayType(WdlIntegerType)), List(
        WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(0), WdlInteger(1))),
        WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(2), WdlInteger(3))))),
      "map" -> WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("foo") -> WdlString("foo"),
        WdlString("bar") -> WdlString("bar"))
      ),
      "map2" -> WdlMap(WdlMapType(WdlStringType, WdlMapType(WdlStringType, WdlStringType)), Map(
        WdlString("foo") ->
          WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlString("foo2") -> WdlString("foo"))),
        WdlString("bar") ->
          WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlString("bar2") -> WdlString("bar")))
      )),
      "map3" -> WdlMap(WdlMapType(WdlStringType, WdlArrayType(WdlIntegerType)), Map(
        WdlString("foo") -> WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(0), WdlInteger(1))),
        WdlString("bar") -> WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(2), WdlInteger(3))))
      ),
      "map4" -> WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("foo[1]") -> WdlString("foo"),
        WdlString("bar[[") -> WdlString("bar"),
        WdlString("baz:qux") -> WdlString("baz:qux")
      ))
    )

    val taskOutputs = wdlValues map { case (k, wv) => TaskOutput(k, wv.wdlType, IgnoredExpression) }

    import WdlValueSimpleton._
    val actual = WdlValueBuilder.toWdlValues(taskOutputs, wdlValues.simplify)
    actual shouldEqual wdlValues
  }
}
