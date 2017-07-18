package wdl4s.wdl.types

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import wdl4s.wdl.WdlNamespaceWithWorkflow
import wdl4s.wdl.values.{WdlArray, WdlInteger, WdlMap, WdlPair, WdlString}

import scala.util.{Failure, Success}

class WdlPairTypeSpec extends FlatSpec with Matchers {

  behavior of "WdlPairType"

  val simplePair = WdlPair(WdlString("a"), WdlInteger(1))

  val stringIntMap = WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
    WdlString("a") -> WdlInteger(1),
    WdlString("b") -> WdlInteger(2),
    WdlString("c") -> WdlInteger(3)
  ))

  val arrayOfPairs = WdlArray(WdlArrayType(WdlPairType(WdlStringType, WdlIntegerType)), Seq(
    WdlPair(WdlString("a"),WdlInteger(1)),
    WdlPair(WdlString("b"),WdlInteger(2)),
    WdlPair(WdlString("c"),WdlInteger(3))
  ))

  val arrayOfPairsOfArrays = WdlArray(WdlArrayType(WdlPairType(WdlArrayType(WdlStringType), WdlArrayType(WdlIntegerType))),
    Seq(
    WdlPair(WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("a"), WdlString("b"))),
      WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(1), WdlInteger(11)))),
    WdlPair(WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c"), WdlString("d"))),
      WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(2), WdlInteger(21)))),
    WdlPair(WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("e"), WdlString("f"))),
      WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(3), WdlInteger(31))))
    )
  )

  val coerceables = Table(
    ("fromValue", "toType", "coercedValue"),
    (
      WdlPair(WdlInteger(1), WdlInteger(2)),
      WdlPairType(WdlStringType, WdlStringType),
      Some(WdlPair(WdlString("1"), WdlString("2")))
    ),

    (
      WdlPair(WdlMap(WdlMapType(WdlIntegerType, WdlStringType), Map(
      WdlInteger(1) -> WdlString("100"),
      WdlInteger(2) -> WdlString("200")
      )), WdlString("300")),
      WdlPairType(WdlMapType(WdlStringType, WdlIntegerType), WdlIntegerType),
      Some(WdlPair(WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
        WdlString("1") -> WdlInteger(100),
        WdlString("2") -> WdlInteger(200)
      )), WdlInteger(300)))
    ),

    (
      WdlPair(WdlMap(WdlMapType(WdlIntegerType, WdlStringType), Map(
        WdlInteger(1) -> WdlString("100"),
        WdlInteger(2) -> WdlString("200")
      )), WdlString("300")),
      WdlPairType(WdlArrayType(WdlStringType), WdlStringType),
      None)

  )

  coerceables foreach { case (fromValue, toType, coercedValue) =>

    val notString = coercedValue map { _ => "" } getOrElse "not "
    val coercionDefined = coercedValue.isDefined

    it should s"${notString}allow coercion from ${fromValue.wdlType.toWdlString} to ${toType.toWdlString}" in {
      toType.isCoerceableFrom(fromValue.wdlType) should be(coercionDefined)
    }

    it should s"generate the expected result when converting from ${fromValue.toWdlString} to ${toType.toWdlString}" in {
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

    WdlArrayType(WdlPairType(WdlStringType, WdlIntegerType)).coerceRawValue(jsArray) match {
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

    WdlArrayType(WdlPairType(WdlArrayType(WdlStringType), WdlArrayType(WdlIntegerType))).coerceRawValue(complexJsArray) match {
      case Success(array) => array shouldEqual arrayOfPairsOfArrays
      case Failure(f) => fail(s"exception while coercing JsObject to WdlPair: $f")
    }
  }

  it should "detect invalid pair construction if missing right or left" in {
    val invalidPair =
      """
        |{ "Left": "a" }
      """.stripMargin.parseJson

    val results = WdlPairType(WdlStringType, WdlIntegerType).coerceRawValue(invalidPair)
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

    val results = WdlPairType(WdlStringType, WdlIntegerType).coerceRawValue(invalidPair)
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

    val results = WdlPairType(WdlStringType, WdlIntegerType).coerceRawValue(invalidPair)
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

    WdlPairType(WdlStringType, WdlIntegerType).coerceRawValue(validCaps) match {
      case Success(array) => array shouldEqual simplePair
      case Failure(f) => fail(s"exception while coercing JsObject to WdlPair: $f")
    }
  }

  it should "return an error if pair types can't be coerced" in {
    val results = WdlPairType(WdlIntegerType, WdlIntegerType).coerceRawValue(JsObject("Left" -> JsString("a")))
    results match {
      case Failure(ex) =>
        ex.getMessage should (startWith("No coercion defined from") and endWith("to 'Pair[Int, Int]'."))
      case Success(_) => fail("Unexpected successful coercion to WdlPair")
    }
  }
}
