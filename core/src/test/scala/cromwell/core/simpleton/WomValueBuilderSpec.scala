package cromwell.core.simpleton

import cromwell.core.simpleton.WomValueBuilderSpec._
import cromwell.util.WomMocks
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wom.callable.Callable.OutputDefinition
import wom.expression.PlaceholderWomExpression
import wom.types.{WomArrayType, WomIntegerType, WomMapType, WomStringType}
import wom.values._

object WomValueBuilderSpec {
  // WdlValueBuilder doesn't care about this expression, but something needs to be passed to the TaskOutput constructor.
  val IgnoredExpression = PlaceholderWomExpression(Set.empty, WomStringType)
}

class WomValueBuilderSpec extends FlatSpec with Matchers with Mockito {

  case class SimpletonConversion(name: String, womValue: WomValue, simpletons: Seq[WomValueSimpleton])
  val simpletonConversions = List(
    SimpletonConversion("foo", WomString("none"), List(WomValueSimpleton("foo", WomString("none")))),
    SimpletonConversion("bar", WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2))), List(WomValueSimpleton("bar[0]", WomInteger(1)), WomValueSimpleton("bar[1]", WomInteger(2)))),
    SimpletonConversion(
      "baz",
      WomArray(WomArrayType(WomArrayType(WomIntegerType)), List(
        WomArray(WomArrayType(WomIntegerType), List(WomInteger(0), WomInteger(1))),
        WomArray(WomArrayType(WomIntegerType), List(WomInteger(2), WomInteger(3))))),
      List(WomValueSimpleton("baz[0][0]", WomInteger(0)), WomValueSimpleton("baz[0][1]", WomInteger(1)), WomValueSimpleton("baz[1][0]", WomInteger(2)), WomValueSimpleton("baz[1][1]", WomInteger(3)))
    ),
    SimpletonConversion(
      "map",
      WomMap(WomMapType(WomStringType, WomStringType), Map(
        WomString("foo") -> WomString("foo"),
        WomString("bar") -> WomString("bar"))),
      List(WomValueSimpleton("map:foo", WomString("foo")), WomValueSimpleton("map:bar", WomString("bar")))
    ),
    SimpletonConversion(
      "mapOfMaps",
      WomMap(WomMapType(WomStringType, WomMapType(WomStringType, WomStringType)), Map(
        WomString("foo") -> WomMap(WomMapType(WomStringType, WomStringType), Map(WomString("foo2") -> WomString("foo"))),
        WomString("bar") ->WomMap(WomMapType(WomStringType, WomStringType), Map(WomString("bar2") -> WomString("bar"))))),
      List(WomValueSimpleton("mapOfMaps:foo:foo2", WomString("foo")), WomValueSimpleton("mapOfMaps:bar:bar2", WomString("bar")))
    ),
    SimpletonConversion(
      "simplePair1",
      WomPair(WomInteger(1), WomString("hello")),
      List(WomValueSimpleton("simplePair1:left", WomInteger(1)), WomValueSimpleton("simplePair1:right", WomString("hello")))
    ),
    SimpletonConversion(
      "simplePair2",
      WomPair(WomString("left"), WomInteger(5)),
      List(WomValueSimpleton("simplePair2:left", WomString("left")), WomValueSimpleton("simplePair2:right", WomInteger(5)))
    ),
    SimpletonConversion(
      "pairOfPairs",
      WomPair(
        WomPair(WomInteger(1), WomString("one")),
        WomPair(WomString("two"), WomInteger(2))),
      List(
        WomValueSimpleton("pairOfPairs:left:left", WomInteger(1)),
        WomValueSimpleton("pairOfPairs:left:right", WomString("one")),
        WomValueSimpleton("pairOfPairs:right:left", WomString("two")),
        WomValueSimpleton("pairOfPairs:right:right", WomInteger(2)))
    ),
    SimpletonConversion(
      "pairOfArrayAndMap",
      WomPair(
        WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2))),
        WomMap(WomMapType(WomStringType, WomIntegerType), Map(WomString("left") -> WomInteger(100), WomString("right") -> WomInteger(200)))),
      List(
        WomValueSimpleton("pairOfArrayAndMap:left[0]", WomInteger(1)),
        WomValueSimpleton("pairOfArrayAndMap:left[1]", WomInteger(2)),
        WomValueSimpleton("pairOfArrayAndMap:right:left", WomInteger(100)),
        WomValueSimpleton("pairOfArrayAndMap:right:right", WomInteger(200)))
    ),
    SimpletonConversion(
      "mapOfArrays",
      WomMap(WomMapType(WomStringType, WomArrayType(WomIntegerType)), Map(
        WomString("foo") -> WomArray(WomArrayType(WomIntegerType), List(WomInteger(0), WomInteger(1))),
        WomString("bar") -> WomArray(WomArrayType(WomIntegerType), List(WomInteger(2), WomInteger(3))))),
      List(WomValueSimpleton("mapOfArrays:foo[0]", WomInteger(0)), WomValueSimpleton("mapOfArrays:foo[1]", WomInteger(1)),
        WomValueSimpleton("mapOfArrays:bar[0]", WomInteger(2)), WomValueSimpleton("mapOfArrays:bar[1]", WomInteger(3)))
    ),
    SimpletonConversion(
      "escapology",
      WomMap(WomMapType(WomStringType, WomStringType), Map(
        WomString("foo[1]") -> WomString("foo"),
        WomString("bar[[") -> WomString("bar"),
        WomString("baz:qux") -> WomString("baz:qux"))),
      List(WomValueSimpleton("escapology:foo\\[1\\]", WomString("foo")),
        WomValueSimpleton("escapology:bar\\[\\[", WomString("bar")),
        WomValueSimpleton("escapology:baz\\:qux", WomString("baz:qux")))
    )
  )

  behavior of "WomValueSimpleton and WdlValueBuilder"

  simpletonConversions foreach { case SimpletonConversion(name, womValue, simpletons) =>
    it should s"decompose WdlValues into simpletons ($name)" in {
      import WomValueSimpleton._

      val map = Map(WomMocks.mockOutputPort(name) -> womValue)
      map.simplify should contain theSameElementsAs simpletons
    }

    it should s"build simpletons back into WdlValues ($name)" in {
      // The task output is used to tell us the type of output we're expecting:
      val outputPort = WomMocks.mockOutputPort(OutputDefinition(name, womValue.womType, IgnoredExpression))
      val taskOutputPorts = List(outputPort)
      val rebuiltValues = WomValueBuilder.toWdlValues(taskOutputPorts, simpletons)
      rebuiltValues.size should be(1)
      rebuiltValues(outputPort) should be(womValue)
    }

  }


  it should "round trip everything together with no losses" in {

    val wdlValues = (simpletonConversions map { case SimpletonConversion(name, womValue, _) => WomMocks.mockOutputPort(name, womValue.womType) -> womValue }).toMap
    val allSimpletons = simpletonConversions flatMap { case SimpletonConversion(_, _, simpletons) => simpletons }

    import WomValueSimpleton._

    val actualSimpletons = wdlValues.simplify
    actualSimpletons should contain theSameElementsAs allSimpletons

    val actual = WomValueBuilder.toWdlValues(wdlValues.keys.toSeq, actualSimpletons)
    actual shouldEqual wdlValues
  }
}
