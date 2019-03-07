package cromwell.core.simpleton

import cromwell.core.simpleton.WomValueBuilderSpec._
import cromwell.core.simpleton.WomValueSimpleton._
import cromwell.util.WomMocks
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wom.callable.Callable.OutputDefinition
import wom.expression.PlaceholderWomExpression
import wom.types.{WomArrayType, WomIntegerType, WomMapType, WomStringType}
import wom.values._

import scala.util.Success

object WomValueBuilderSpec {
  // WdlValueBuilder doesn't care about this expression, but something needs to be passed to the TaskOutput constructor.
  val IgnoredExpression = PlaceholderWomExpression(Set.empty, WomStringType)
}

class WomValueBuilderSpec extends FlatSpec with Matchers with Mockito {

  case class SimpletonConversion(name: String, womValue: WomValue, simpletons: Seq[WomValueSimpleton])
  val simpletonConversions = List(
    SimpletonConversion("foo", WomString("none"), List(WomValueSimpleton("foo", WomString("none")))),
    SimpletonConversion("bar", WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2))), List(WomValueSimpleton("bar[0]", WomInteger(1)), WomValueSimpleton("bar[1]", WomInteger(2)))),
    SimpletonConversion("empty_array", WomArray(WomArrayType(WomIntegerType), List.empty), List()),
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
    ),
    SimpletonConversion(
      "flat_object",
      WomObject(Map(
        "a" -> WomString("aardvark"),
        "b" -> WomInteger(25),
        "c" -> WomBoolean(false)
      )),
      List(WomValueSimpleton("flat_object:a", WomString("aardvark")),
        WomValueSimpleton("flat_object:b", WomInteger(25)),
        WomValueSimpleton("flat_object:c", WomBoolean(false)))
    ),
    SimpletonConversion(
      "object_with_array",
      WomObject(Map(
        "a" -> WomArray(WomArrayType(WomStringType), Seq(WomString("aardvark"), WomString("beetle")))
      )),
      List(WomValueSimpleton("object_with_array:a[0]", WomString("aardvark")),
        WomValueSimpleton("object_with_array:a[1]", WomString("beetle")))
    ),
    SimpletonConversion(
      "object_with_object",
      WomObject(Map(
        "a" -> WomObject(Map(
          "aa" -> WomArray(WomArrayType(WomStringType), Seq(WomString("aardvark"), WomString("aaron"))),
          "ab" -> WomArray(WomArrayType(WomStringType), Seq(WomString("abacus"), WomString("a bee")))
        )),
        "b" -> WomObject(Map(
          "ba" -> WomArray(WomArrayType(WomStringType), Seq(WomString("baa"), WomString("battle"))),
          "bb" -> WomArray(WomArrayType(WomStringType), Seq(WomString("bbrrrr"), WomString("bb gun")))
        ))
      )),
      List(
        WomValueSimpleton("object_with_object:a:aa[0]", WomString("aardvark")),
        WomValueSimpleton("object_with_object:a:aa[1]", WomString("aaron")),
        WomValueSimpleton("object_with_object:a:ab[0]", WomString("abacus")),
        WomValueSimpleton("object_with_object:a:ab[1]", WomString("a bee")),
        WomValueSimpleton("object_with_object:b:ba[0]", WomString("baa")),
        WomValueSimpleton("object_with_object:b:ba[1]", WomString("battle")),
        WomValueSimpleton("object_with_object:b:bb[0]", WomString("bbrrrr")),
        WomValueSimpleton("object_with_object:b:bb[1]", WomString("bb gun")),
      )
    ),
    /*
      * Wom object representing a directory listing
      *  - a single file
      *  - a "maybe populated file" with some properties (checksum etc..) and secondary files:
      *         - another single file
      *         - a directory listing a single file
      *         - an unlisted directory
      *         - a glob file
      *  - an unlisted directory
      *  - a glob file
      *  
      *  Note: glob files technically are never simpletonized but as WomFiles they *can* be
     */
    SimpletonConversion(
      "directory",
      WomMaybeListedDirectory(
        Option("outerValueName"), 
        Option(List(
          WomSingleFile("outerSingleFile"),
          WomMaybeListedDirectory(Option("innerValueName"), Option(List(WomSingleFile("innerSingleFile")))),
          WomMaybePopulatedFile(
            Option("populatedInnerValueName"),
            Option("innerChecksum"),
            Option(10L),
            Option("innerFormat"),
            Option("innerContents"),
            List(
              WomSingleFile("populatedInnerSingleFile"),
              WomMaybeListedDirectory(Option("innerDirectoryValueName"), Option(List(WomSingleFile("innerDirectorySingleFile")))),
              WomUnlistedDirectory("innerUnlistedDirectory"),
              WomGlobFile("innerGlobFile")
            )
          ),
          WomUnlistedDirectory("outerUnlistedDirectory"),
          WomGlobFile("outerGlobFile")
      ))),
      List(
        WomValueSimpleton("directory:class", WomString("Directory")),
        WomValueSimpleton("directory:value", WomString("outerValueName")),
        
        WomValueSimpleton("directory:listing[0]", WomSingleFile("outerSingleFile")),
        
        WomValueSimpleton("directory:listing[1]:class", WomString("Directory")),
        WomValueSimpleton("directory:listing[1]:value", WomString("innerValueName")),
        WomValueSimpleton("directory:listing[1]:listing[0]", WomSingleFile("innerSingleFile")),
        
        WomValueSimpleton("directory:listing[2]:class", WomString("File")),
        WomValueSimpleton("directory:listing[2]:value", WomString("populatedInnerValueName")),
        
        WomValueSimpleton("directory:listing[2]:checksum", WomString("innerChecksum")),
        WomValueSimpleton("directory:listing[2]:size", WomInteger(10)),
        WomValueSimpleton("directory:listing[2]:format", WomString("innerFormat")),
        WomValueSimpleton("directory:listing[2]:contents", WomString("innerContents")),
        WomValueSimpleton("directory:listing[2]:secondaryFiles[0]", WomSingleFile("populatedInnerSingleFile")),
        WomValueSimpleton("directory:listing[2]:secondaryFiles[1]:class", WomString("Directory")),
        WomValueSimpleton("directory:listing[2]:secondaryFiles[1]:value", WomString("innerDirectoryValueName")),
        WomValueSimpleton("directory:listing[2]:secondaryFiles[1]:listing[0]", WomSingleFile("innerDirectorySingleFile")),
        WomValueSimpleton("directory:listing[2]:secondaryFiles[2]", WomUnlistedDirectory("innerUnlistedDirectory")),
        WomValueSimpleton("directory:listing[2]:secondaryFiles[3]", WomGlobFile("innerGlobFile")),
        
        WomValueSimpleton("directory:listing[3]", WomUnlistedDirectory("outerUnlistedDirectory")),
        WomValueSimpleton("directory:listing[4]", WomGlobFile("outerGlobFile"))
      )
    )
  )

  behavior of "WomValueSimpleton and WdlValueBuilder"

  simpletonConversions foreach { case SimpletonConversion(name, womValue, expectedSimpletons) =>
    it should s"decompose WdlValues into simpletons ($name)" in {

      val map = Map(WomMocks.mockOutputPort(name) -> womValue)
      assertSimpletonsEqual(expectedSimpletons, map.simplify)
    }

    it should s"build simpletons back into WdlValues ($name)" in {
      // The task output is used to tell us the type of output we're expecting:
      val outputPort = WomMocks.mockOutputPort(OutputDefinition(name, womValue.womType, IgnoredExpression))
      val taskOutputPorts = List(outputPort)
      val rebuiltValues = WomValueBuilder.toWomValues(taskOutputPorts, expectedSimpletons)
      rebuiltValues.size should be(1)
      rebuiltValues(outputPort) should be(womValue)
    }
  }

  it should "round trip everything together with no losses" in {

    val wdlValues = (simpletonConversions map { case SimpletonConversion(name, womValue, _) => WomMocks.mockOutputPort(name, womValue.womType) -> womValue }).toMap
    val allSimpletons = simpletonConversions flatMap { case SimpletonConversion(_, _, simpletons) => simpletons }

    val actualSimpletons = wdlValues.simplify
    assertSimpletonsEqual(allSimpletons, actualSimpletons)

    val actual = WomValueBuilder.toWomValues(wdlValues.keys.toSeq, actualSimpletons)
    actual shouldEqual wdlValues
  }

  // We won't get exactly the same thing back when we reconstruct a Map from inside an object, but it should be
  // coerceable back into the original type:
  it should "decompose then reconstruct a map in an object into a coerceable value" in {

    val aMap = WomMap(WomMapType(WomStringType, WomArrayType(WomStringType)), Map(
      WomString("aa") -> WomArray(WomArrayType(WomStringType), Seq(WomString("aardvark"), WomString("aaron"))),
      WomString("ab") -> WomArray(WomArrayType(WomStringType), Seq(WomString("abacus"), WomString("a bee")))
    ))

    val bMap = WomMap(WomMapType(WomStringType, WomArrayType(WomStringType)), Map(
      WomString("ba") -> WomArray(WomArrayType(WomStringType), Seq(WomString("baa"), WomString("battle"))),
      WomString("bb") -> WomArray(WomArrayType(WomStringType), Seq(WomString("bbrrrr"), WomString("bb gun")))
    ))

    val initial = WomObject(Map("a" -> aMap, "b" -> bMap ))

    val map = Map(WomMocks.mockOutputPort("map_in_object") -> initial)

    val actualSimpletons = map.simplify
    assertSimpletonsEqual(
      List(
        WomValueSimpleton("map_in_object:a:aa[0]", WomString("aardvark")),
        WomValueSimpleton("map_in_object:a:aa[1]", WomString("aaron")),
        WomValueSimpleton("map_in_object:a:ab[0]", WomString("abacus")),
        WomValueSimpleton("map_in_object:a:ab[1]", WomString("a bee")),
        WomValueSimpleton("map_in_object:b:ba[0]", WomString("baa")),
        WomValueSimpleton("map_in_object:b:ba[1]", WomString("battle")),
        WomValueSimpleton("map_in_object:b:bb[0]", WomString("bbrrrr")),
        WomValueSimpleton("map_in_object:b:bb[1]", WomString("bb gun")),
      ),
      actualSimpletons)

    // Reconstruct:
    val outputPort = WomMocks.mockOutputPort(OutputDefinition("map_in_object", initial.womType, IgnoredExpression))
    val taskOutputPorts = List(outputPort)
    val rebuiltValues = WomValueBuilder.toWomValues(taskOutputPorts, actualSimpletons)

    rebuiltValues.size should be(1)
    val rebuiltObject = rebuiltValues.head._2
    rebuiltObject match {
      case o: WomObject =>
        aMap.womType.coerceRawValue(o.values("a")) should be(Success(aMap))
        bMap.womType.coerceRawValue(o.values("b")) should be(Success(bMap))
      case other => fail(s"Expected reconstruction to Object but got ${other.womType.stableName}")
    }
  }

  def assertSimpletonsEqual(expectedSimpletons: Iterable[WomValueSimpleton], actualSimpletons: Iterable[WomValueSimpleton]) = {

    // Sanity check, make sure we don't lose anything when we "toSet":
    actualSimpletons.toSet should contain theSameElementsAs actualSimpletons
    actualSimpletons.toSet.size should be(actualSimpletons.size)
    expectedSimpletons.toSet should contain theSameElementsAs expectedSimpletons
    expectedSimpletons.toSet.size should be(expectedSimpletons.size)

    if (actualSimpletons.toSet != expectedSimpletons.toSet) {
      val unexpecteds = actualSimpletons.toSet.diff(expectedSimpletons.toSet)
      val unfounds = expectedSimpletons.toSet.diff(actualSimpletons.toSet)
      fail(
        s"""Actual simpletons did not meet expectations
           |Total found / total expected: ${actualSimpletons.size} / ${expectedSimpletons.size}
           |Found but not expected: ${unexpecteds.mkString(", ")}
           |Expected but not found: ${unfounds.mkString(", ")}
           |""".stripMargin
      )
    }
  }
}
