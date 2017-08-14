package wdl4s.wdl.types

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.values.{WdlArray, WdlFile, WdlInteger, WdlMap, WdlOptionalValue, WdlString}

import scala.util.{Failure, Success}

class WdlOptionalTypeSpec extends FlatSpec with Matchers {

  import TableDrivenPropertyChecks._

  behavior of "WdlOptionalType"

  val innerCoercionDefinedFromTable = Table(
    ("from type", "to type", "coerceableFrom"),
    ( WdlOptionalType(WdlStringType), WdlOptionalType(WdlFileType), true ),
    ( WdlOptionalType(WdlOptionalType(WdlStringType)), WdlOptionalType(WdlOptionalType(WdlFileType)), true ),
    ( WdlOptionalType(WdlIntegerType), WdlOptionalType(WdlBooleanType), false ),
    ( WdlOptionalType(WdlOptionalType(WdlIntegerType)), WdlOptionalType(WdlOptionalType(WdlBooleanType)), false )
  )

  val boxingCoercionDefinedFromTable = Table(
    ("from type", "to type", "coerceableFrom"),
    ( WdlStringType, WdlOptionalType(WdlStringType), true ),
    ( WdlStringType, WdlOptionalType(WdlOptionalType(WdlStringType)), true ),
    ( WdlStringType, WdlOptionalType(WdlFileType), true ),
    ( WdlOptionalType(WdlStringType), WdlOptionalType(WdlOptionalType(WdlStringType)), true ),
    ( WdlOptionalType(WdlStringType), WdlStringType, false )
  )

  innerCoercionDefinedFromTable foreach { case (from, to, expectation) =>
    val notString = if (!expectation) "not" else ""
    it should s"$notString allow coercion from ${from.toWdlString} to ${to.toWdlString}" in {
      to.isCoerceableFrom(from) should be(expectation)
    }
  }

  boxingCoercionDefinedFromTable foreach { case (from, to, expectation) =>
    val notString = if (!expectation) "not" else ""
    it should s"$notString allow boxing coercion from ${from.toWdlString} to ${to.toWdlString}" in {
      to.isCoerceableFrom(from) should be(expectation)
    }
  }

  val innerCoercionTable = Table(
    ("from value", "to type", "expectation"),
    ( WdlOptionalValue(WdlString("3")), WdlOptionalType(WdlIntegerType), WdlOptionalValue(WdlInteger(3)) ),
    ( WdlOptionalValue(WdlOptionalValue(WdlString("3"))), WdlOptionalType(WdlOptionalType(WdlIntegerType)), WdlOptionalValue(WdlOptionalValue(WdlInteger(3))) )
  )

  innerCoercionTable foreach { case (from, to, expectation) =>
    it should s"coerce correctly from ${from.toWdlString} to ${to.toWdlString}" in {
      to.coerceRawValue(from) should be(Success(expectation))
    }
  }

  val boxingCoercionTable = Table(
    ("from value", "to type", "expectation"),
    ( WdlString("boxme1"), WdlOptionalType(WdlStringType), WdlOptionalValue(WdlString("boxme1")) ),
    ( WdlString("boxme2"), WdlOptionalType(WdlOptionalType(WdlStringType)), WdlOptionalValue(WdlOptionalValue(WdlString("boxme2"))) ),
    ( WdlString("boxme3"), WdlOptionalType(WdlFileType), WdlOptionalValue(WdlFile("boxme3")) ),
    ( WdlOptionalValue(WdlString("boxme4")), WdlOptionalType(WdlOptionalType(WdlStringType)), WdlOptionalValue(WdlOptionalValue(WdlString("boxme4"))) )
  )

  boxingCoercionTable foreach { case (from, to, expectation) =>
    it should s"coerce correctly from ${from.toWdlString} to ${to.toWdlString}" in {
      to.coerceRawValue(from) match {
        case Success(coerced) => coerced should be(expectation)
        case Failure(t) => fail("Unexpected coercion failure: " + t)
      }
    }
  }

  it should "coerce JsNull to empty value" in {
    import spray.json._

    // Literal null
    WdlOptionalType(WdlIntegerType).coerceRawValue(JsNull) shouldBe Success(WdlOptionalValue.none(WdlIntegerType))
    // Null in array
    val jsArray = "[1, 2, null, 4]".parseJson
    val arrayType = WdlArrayType(WdlOptionalType(WdlIntegerType))
    arrayType.coerceRawValue(jsArray) shouldBe
      Success(WdlArray(arrayType,
        List(WdlOptionalValue(WdlInteger(1)), WdlOptionalValue(WdlInteger(2)), WdlOptionalValue.none(WdlIntegerType), WdlOptionalValue(WdlInteger(4)))
      ))
    // Null in object
    val jsObject =
      """
        |{
        |  "one": 1,
        |  "two": null,
        |  "three": 3
        |}
      """.stripMargin.parseJson
    val mapType = WdlMapType(WdlStringType, WdlOptionalType(WdlIntegerType))
    mapType.coerceRawValue(jsObject) shouldBe
      Success(WdlMap(mapType,
        Map(
          WdlString("one") -> WdlOptionalValue(WdlInteger(1)),
          WdlString("two") -> WdlOptionalValue.none(WdlIntegerType),
          WdlString("three") ->WdlOptionalValue(WdlInteger(3))
        )
      ))
  }

}
