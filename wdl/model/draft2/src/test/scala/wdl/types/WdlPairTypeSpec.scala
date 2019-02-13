package wdl.types

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import wdl.draft2.model.WdlNamespaceWithWorkflow
import wom.types._
import wom.values.{WomInteger, WomMap, WomPair, WomString, _}

import scala.util.{Failure, Success}

class WdlPairTypeSpec extends FlatSpec with Matchers {

  behavior of "WdlPairType"

  val simplePair = WomPair(WomString("a"), WomInteger(1))

  val stringIntMap = WomMap(WomMapType(WomStringType, WomIntegerType), Map(
    WomString("a") -> WomInteger(1),
    WomString("b") -> WomInteger(2),
    WomString("c") -> WomInteger(3)
  ))

  val arrayOfPairs = WomArray(WomArrayType(WomPairType(WomStringType, WomIntegerType)), Seq(
    WomPair(WomString("a"),WomInteger(1)),
    WomPair(WomString("b"),WomInteger(2)),
    WomPair(WomString("c"),WomInteger(3))
  ))

  val arrayOfPairsOfArrays = WomArray(WomArrayType(WomPairType(WomArrayType(WomStringType), WomArrayType(WomIntegerType))),
    Seq(
    WomPair(WomArray(WomArrayType(WomStringType), Seq(WomString("a"), WomString("b"))),
      WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(1), WomInteger(11)))),
    WomPair(WomArray(WomArrayType(WomStringType), Seq(WomString("c"), WomString("d"))),
      WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(2), WomInteger(21)))),
    WomPair(WomArray(WomArrayType(WomStringType), Seq(WomString("e"), WomString("f"))),
      WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(3), WomInteger(31))))
    )
  )

  val coerceables = Table(
    ("fromValue", "toType", "coercedValue"),
    (
      WomPair(WomInteger(1), WomInteger(2)),
      WomPairType(WomStringType, WomStringType),
      Some(WomPair(WomString("1"), WomString("2")))
    ),

    (
      WomPair(WomMap(WomMapType(WomIntegerType, WomStringType), Map(
      WomInteger(1) -> WomString("100"),
      WomInteger(2) -> WomString("200")
      )), WomString("300")),
      WomPairType(WomMapType(WomStringType, WomIntegerType), WomIntegerType),
      Some(WomPair(WomMap(WomMapType(WomStringType, WomIntegerType), Map(
        WomString("1") -> WomInteger(100),
        WomString("2") -> WomInteger(200)
      )), WomInteger(300)))
    ),

    (
      WomPair(WomMap(WomMapType(WomIntegerType, WomStringType), Map(
        WomInteger(1) -> WomString("100"),
        WomInteger(2) -> WomString("200")
      )), WomString("300")),
      WomPairType(WomArrayType(WomStringType), WomStringType),
      None)

  )

  coerceables foreach { case (fromValue, toType, coercedValue) =>

    val notString = coercedValue map { _ => "" } getOrElse "not "
    val coercionDefined = coercedValue.isDefined

    it should s"${notString}allow coercion from ${fromValue.womType.stableName} to ${toType.stableName}" in {
      toType.isCoerceableFrom(fromValue.womType) should be(coercionDefined)
    }

    it should s"generate the expected result when converting from ${fromValue.toWomString} to ${toType.stableName}" in {
      val actual = toType.coerceRawValue(fromValue)
      (actual, coercedValue) match {
        case (Success(actualValue), Some(expectedValue)) => actualValue should be(expectedValue)
        case (Success(actualValue), None) => fail(s"Coercion should have failed but instead got $actualValue")
        case (Failure(_), None) => // Correctly failed to coerce
        case (Failure(t), Some(expectedValue)) => fail(s"Expected coercion to produce $expectedValue but instead got exception $t")
      }
    }
  }

  it should "Allow use of pairs in call input mappings" in {
    val wdl =
      """
        |task t {
        |  File f1
        |  File f2
        |  File f3
        |  File f4
        |  File f5
        |  command {...}
        |  output {
        |    File o1 = f1
        |    File o2 = f2
        |    File o3 = f3
        |    File o4 = f4
        |    File o5 = f5
        |  }
        |}
        |
        |workflow w {
        |  Array[Pair[File, Array[File]]] arr = [("left_file", ["right_file1", "right_file2"])]
        |  Pair[String, String] p = ("hello", "world")
        |
        |  scatter (i in arr) {
        |    call t {
        |      input:
        |       f1 = i.left,
        |       f2 = i.right[0],
        |       f3 = i.right[1],
        |       f4 = p.left,
        |       f5 = p.right
        |    }
        |  }
        |}
      """.stripMargin

    noException should be thrownBy WdlNamespaceWithWorkflow.load(wdl, Seq.empty)
  }

  it should "coerce a JsArray into a WdlArray of WdlPairs" in {
    val jsArray =
      """
        |[ { "Left": "a", "Right": 1 },
        |  { "Left": "b", "Right": 2 },
        |  { "Left": "c", "Right": 3 }
        |]
      """.stripMargin.parseJson

    WomArrayType(WomPairType(WomStringType, WomIntegerType)).coerceRawValue(jsArray) match {
      case Success(array) => array shouldEqual arrayOfPairs
      case Failure(f) => fail(s"exception while coercing JsObject: $f")
    }
  }

  it should "coerce a complex JsArray into a WdlArray of WdlPairs of WdlArrays" in {
    val complexJsArray =
      """
        |[ { "Left": ["a", "b"], "Right": [1, 11] },
        |  { "Left": ["c", "d"], "Right": [2, 21] },
        |  { "Left": ["e", "f"], "Right": [3, 31] }
        |]
      """.stripMargin.parseJson

    WomArrayType(WomPairType(WomArrayType(WomStringType), WomArrayType(WomIntegerType))).coerceRawValue(complexJsArray) match {
      case Success(array) => array shouldEqual arrayOfPairsOfArrays
      case Failure(f) => fail(s"exception while coercing JsObject to WdlPair: $f")
    }
  }

  it should "detect invalid pair construction if missing right or left" in {
    val invalidPair =
      """
        |{ "Left": "a" }
      """.stripMargin.parseJson

    val results = WomPairType(WomStringType, WomIntegerType).coerceRawValue(invalidPair)
    results match {
      case Failure(ex) =>
        ex.getMessage should (startWith("No coercion defined from") and endWith("to 'Pair[String, Int]'."))
      case Success(_) => fail("Unexpected successful coercion to WdlPair")
    }
  }

  it should "detect invalid pair construction if JsObject is of size >2" in {
    val invalidPair =
      """
        |{ "Left": "a", "Middle": "notAllowed", "Right": 1 }
      """.stripMargin.parseJson

    val results = WomPairType(WomStringType, WomIntegerType).coerceRawValue(invalidPair)
    results match {
      case Failure(ex) =>
        ex.getMessage should (startWith("No coercion defined from") and endWith("to 'Pair[String, Int]'."))
      case Success(_) => fail("Unexpected successful coercion to WdlPair")
    }
  }

  it should "detect invalid pair construction if Left and/or Right are undefined in the JsObject" in {
    val invalidPair =
      """
        |{ "notLeft": "a", "notRight": 1 }
      """.stripMargin.parseJson

    val results = WomPairType(WomStringType, WomIntegerType).coerceRawValue(invalidPair)
    results match {
      case Failure(ex) =>
        ex.getMessage should (startWith("Failed to coerce") and endWith("requires for Right/Left value(s) to be defined.)"))
      case Success(_) => fail("Unexpected successful coercion to WdlPair")
    }
  }

  it should "coerce a JsObject to a WdlPair regardless of canonical capitalization" in {
    val validCaps =
      """
        |{ "LefT": "a", "right": 1 }
      """.stripMargin.parseJson

    WomPairType(WomStringType, WomIntegerType).coerceRawValue(validCaps) match {
      case Success(array) => array shouldEqual simplePair
      case Failure(f) => fail(s"exception while coercing JsObject to WdlPair: $f")
    }
  }

  it should "return an error if pair types can't be coerced" in {
    val results = WomPairType(WomIntegerType, WomIntegerType).coerceRawValue(JsObject("Left" -> JsString("a")))
    results match {
      case Failure(ex) =>
        ex.getMessage should (startWith("No coercion defined from") and endWith("to 'Pair[Int, Int]'."))
      case Success(_) => fail("Unexpected successful coercion to WdlPair")
    }
  }
}
