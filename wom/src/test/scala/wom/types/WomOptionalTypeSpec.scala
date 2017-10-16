package wom.types

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.values._

import scala.util.{Failure, Success}

class WomOptionalTypeSpec extends FlatSpec with Matchers {

  import TableDrivenPropertyChecks._

  behavior of "WomOptionalType"

  val innerCoercionDefinedFromTable = Table(
    ("from type", "to type", "coerceableFrom"),
    ( WomOptionalType(WomStringType), WomOptionalType(WomFileType), true ),
    ( WomOptionalType(WomOptionalType(WomStringType)), WomOptionalType(WomOptionalType(WomFileType)), true ),
    ( WomOptionalType(WomIntegerType), WomOptionalType(WomBooleanType), false ),
    ( WomOptionalType(WomOptionalType(WomIntegerType)), WomOptionalType(WomOptionalType(WomBooleanType)), false )
  )

  val boxingCoercionDefinedFromTable = Table(
    ("from type", "to type", "coerceableFrom"),
    ( WomStringType, WomOptionalType(WomStringType), true ),
    ( WomStringType, WomOptionalType(WomOptionalType(WomStringType)), true ),
    ( WomStringType, WomOptionalType(WomFileType), true ),
    ( WomOptionalType(WomStringType), WomOptionalType(WomOptionalType(WomStringType)), true ),
    ( WomOptionalType(WomStringType), WomStringType, false )
  )

  innerCoercionDefinedFromTable foreach { case (from, to, expectation) =>
    val notString = if (!expectation) "not" else ""
    it should s"$notString allow coercion from ${from.toDisplayString} to ${to.toDisplayString}" in {
      to.isCoerceableFrom(from) should be(expectation)
    }
  }

  boxingCoercionDefinedFromTable foreach { case (from, to, expectation) =>
    val notString = if (!expectation) "not" else ""
    it should s"$notString allow boxing coercion from ${from.toDisplayString} to ${to.toDisplayString}" in {
      to.isCoerceableFrom(from) should be(expectation)
    }
  }

  val innerCoercionTable = Table(
    ("from value", "to type", "expectation"),
    ( WomOptionalValue(WomString("3")), WomOptionalType(WomIntegerType), WomOptionalValue(WomInteger(3)) ),
    ( WomOptionalValue(WomOptionalValue(WomString("3"))), WomOptionalType(WomOptionalType(WomIntegerType)), WomOptionalValue(WomOptionalValue(WomInteger(3))) )
  )

  innerCoercionTable foreach { case (from, to, expectation) =>
    it should s"coerce correctly from ${from.toWomString} to ${to.toDisplayString}" in {
      to.coerceRawValue(from) should be(Success(expectation))
    }
  }

  val boxingCoercionTable = Table(
    ("from value", "to type", "expectation"),
    ( WomString("boxme1"), WomOptionalType(WomStringType), WomOptionalValue(WomString("boxme1")) ),
    ( WomString("boxme2"), WomOptionalType(WomOptionalType(WomStringType)), WomOptionalValue(WomOptionalValue(WomString("boxme2"))) ),
    ( WomString("boxme3"), WomOptionalType(WomFileType), WomOptionalValue(WomFile("boxme3")) ),
    ( WomOptionalValue(WomString("boxme4")), WomOptionalType(WomOptionalType(WomStringType)), WomOptionalValue(WomOptionalValue(WomString("boxme4"))) )
  )

  boxingCoercionTable foreach { case (from, to, expectation) =>
    it should s"coerce correctly from ${from.toWomString} to ${to.toDisplayString}" in {
      to.coerceRawValue(from) match {
        case Success(coerced) => coerced should be(expectation)
        case Failure(t) => fail("Unexpected coercion failure: " + t)
      }
    }
  }

  it should "coerce JsNull to empty value" in {
    import spray.json._

    // Literal null
    WomOptionalType(WomIntegerType).coerceRawValue(JsNull) shouldBe Success(WomOptionalValue.none(WomIntegerType))
    // Null in array
    val jsArray = "[1, 2, null, 4]".parseJson
    val arrayType = WomArrayType(WomOptionalType(WomIntegerType))
    arrayType.coerceRawValue(jsArray) shouldBe
      Success(WomArray(arrayType,
        List(WomOptionalValue(WomInteger(1)), WomOptionalValue(WomInteger(2)), WomOptionalValue.none(WomIntegerType), WomOptionalValue(WomInteger(4)))
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
    val mapType = WomMapType(WomStringType, WomOptionalType(WomIntegerType))
    mapType.coerceRawValue(jsObject) shouldBe
      Success(WomMap(mapType,
        Map(
          WomString("one") -> WomOptionalValue(WomInteger(1)),
          WomString("two") -> WomOptionalValue.none(WomIntegerType),
          WomString("three") ->WomOptionalValue(WomInteger(3))
        )
      ))
  }

}
