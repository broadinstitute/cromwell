package cromwell.core.simpleton

import cromwell.core.simpleton.WdlValueBuilderSpec._
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s.parser.WdlParser.Ast
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlMapType, WdlStringType}
import wdl4s.values.{WdlArray, WdlInteger, WdlMap, WdlPair, WdlString, WdlValue}
import wdl4s.{TaskOutput, WdlExpression}

object WdlValueBuilderSpec {
  // WdlValueBuilder doesn't care about this expression, but something needs to be passed to the TaskOutput constructor.
  val IgnoredExpression = WdlExpression.fromString(""" "" """)
}

class WdlValueBuilderSpec extends FlatSpec with Matchers with Mockito {

  case class SimpletonConversion(name: String, wdlValue: WdlValue, simpletons: Seq[WdlValueSimpleton])
  val simpletonConversions = List(
    SimpletonConversion("foo", WdlString("none"), List(WdlValueSimpleton("foo", WdlString("none")))),
    SimpletonConversion("bar", WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(2))), List(WdlValueSimpleton("bar[0]", WdlInteger(1)), WdlValueSimpleton("bar[1]", WdlInteger(2)))),
    SimpletonConversion(
      "baz",
      WdlArray(WdlArrayType(WdlArrayType(WdlIntegerType)), List(
        WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(0), WdlInteger(1))),
        WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(2), WdlInteger(3))))),
      List(WdlValueSimpleton("baz[0][0]", WdlInteger(0)), WdlValueSimpleton("baz[0][1]", WdlInteger(1)), WdlValueSimpleton("baz[1][0]", WdlInteger(2)), WdlValueSimpleton("baz[1][1]", WdlInteger(3)))
    ),
    SimpletonConversion(
      "map",
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("foo") -> WdlString("foo"),
        WdlString("bar") -> WdlString("bar"))),
      List(WdlValueSimpleton("map:foo", WdlString("foo")), WdlValueSimpleton("map:bar", WdlString("bar")))
    ),
    SimpletonConversion(
      "mapOfMaps",
      WdlMap(WdlMapType(WdlStringType, WdlMapType(WdlStringType, WdlStringType)), Map(
        WdlString("foo") -> WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlString("foo2") -> WdlString("foo"))),
        WdlString("bar") ->WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlString("bar2") -> WdlString("bar"))))),
      List(WdlValueSimpleton("mapOfMaps:foo:foo2", WdlString("foo")), WdlValueSimpleton("mapOfMaps:bar:bar2", WdlString("bar")))
    ),
    SimpletonConversion(
      "simplePair1",
      WdlPair(WdlInteger(1), WdlString("hello")),
      List(WdlValueSimpleton("simplePair1:left", WdlInteger(1)), WdlValueSimpleton("simplePair1:right", WdlString("hello")))
    ),
    SimpletonConversion(
      "simplePair2",
      WdlPair(WdlString("left"), WdlInteger(5)),
      List(WdlValueSimpleton("simplePair2:left", WdlString("left")), WdlValueSimpleton("simplePair2:right", WdlInteger(5)))
    ),
    SimpletonConversion(
      "pairOfPairs",
      WdlPair(
        WdlPair(WdlInteger(1), WdlString("one")),
        WdlPair(WdlString("two"), WdlInteger(2))),
      List(
        WdlValueSimpleton("pairOfPairs:left:left", WdlInteger(1)),
        WdlValueSimpleton("pairOfPairs:left:right", WdlString("one")),
        WdlValueSimpleton("pairOfPairs:right:left", WdlString("two")),
        WdlValueSimpleton("pairOfPairs:right:right", WdlInteger(2)))
    ),
    SimpletonConversion(
      "pairOfArrayAndMap",
      WdlPair(
        WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(2))),
        WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(WdlString("left") -> WdlInteger(100), WdlString("right") -> WdlInteger(200)))),
      List(
        WdlValueSimpleton("pairOfArrayAndMap:left[0]", WdlInteger(1)),
        WdlValueSimpleton("pairOfArrayAndMap:left[1]", WdlInteger(2)),
        WdlValueSimpleton("pairOfArrayAndMap:right:left", WdlInteger(100)),
        WdlValueSimpleton("pairOfArrayAndMap:right:right", WdlInteger(200)))
    ),
    SimpletonConversion(
      "mapOfArrays",
      WdlMap(WdlMapType(WdlStringType, WdlArrayType(WdlIntegerType)), Map(
        WdlString("foo") -> WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(0), WdlInteger(1))),
        WdlString("bar") -> WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(2), WdlInteger(3))))),
      List(WdlValueSimpleton("mapOfArrays:foo[0]", WdlInteger(0)), WdlValueSimpleton("mapOfArrays:foo[1]", WdlInteger(1)),
        WdlValueSimpleton("mapOfArrays:bar[0]", WdlInteger(2)), WdlValueSimpleton("mapOfArrays:bar[1]", WdlInteger(3)))
    ),
    SimpletonConversion(
      "escapology",
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("foo[1]") -> WdlString("foo"),
        WdlString("bar[[") -> WdlString("bar"),
        WdlString("baz:qux") -> WdlString("baz:qux"))),
      List(WdlValueSimpleton("escapology:foo\\[1\\]", WdlString("foo")),
        WdlValueSimpleton("escapology:bar\\[\\[", WdlString("bar")),
        WdlValueSimpleton("escapology:baz\\:qux", WdlString("baz:qux")))
    )
  )

  behavior of "WdlValueSimpleton and WdlValueBuilder"

  simpletonConversions foreach { case SimpletonConversion(name, wdlValue, simpletons) =>
    it should s"decompose WdlValues into simpletons ($name)" in {
      import WdlValueSimpleton._

      val map = Map(name -> wdlValue)
      map.simplify should contain theSameElementsAs simpletons
    }

    it should s"build simpletons back into WdlValues ($name)" in {
      // The task output is used to tell us the type of output we're expecting:
      val taskOutputs = List(TaskOutput(name, wdlValue.wdlType, IgnoredExpression, mock[Ast], None))
      val rebuiltValues = WdlValueBuilder.toWdlValues(taskOutputs, simpletons)
      rebuiltValues.size should be(1)
      rebuiltValues(name) should be(wdlValue)
    }

  }


  it should "round trip everything together with no losses" in {

    val wdlValues = (simpletonConversions map { case SimpletonConversion(name, wdlValue, simpletons) => name -> wdlValue }).toMap
    val taskOutputs = wdlValues map { case (k, wv) => TaskOutput(k, wv.wdlType, IgnoredExpression, mock[Ast], None) }
    val allSimpletons = simpletonConversions flatMap { case SimpletonConversion(name, wdlValue, simpletons) => simpletons }

    import WdlValueSimpleton._

    val actualSimpletons = wdlValues.simplify
    actualSimpletons should contain theSameElementsAs allSimpletons

    val actual = WdlValueBuilder.toWdlValues(taskOutputs, actualSimpletons)
    actual shouldEqual wdlValues
  }
}
