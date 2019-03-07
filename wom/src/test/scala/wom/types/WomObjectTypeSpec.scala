package wom.types

import org.scalatest.{FlatSpec, Matchers}
import wom.values.{WomInteger, WomMap, WomObject, WomString}

import scala.util.{Failure, Success}

class WomObjectTypeSpec extends FlatSpec with Matchers {
  val abcObject = WomObject(Map(
    "a" -> WomString("one"),
    "b" -> WomString("two"),
    "c" -> WomString("three")
  ))

  val coerceableMap = WomMap(WomMapType(WomStringType, WomStringType), Map(
    WomString("a") -> WomString("one"),
    WomString("b") -> WomString("two"),
    WomString("c") -> WomString("three"))
  )

  val abcMixedTypedObject = WomObject(Map(
    "a" -> WomString("one"),
    "b" -> WomInteger(2),
    "c" -> WomString("three")
  ))

  val coerceableTypedObject = WomObject.withTypeUnsafe(Map(
    "a" -> WomString("one"),
    "b" -> WomInteger(2),
    "c" -> WomString("three")),
    WomCompositeType(Map("a" -> WomStringType, "b" -> WomIntegerType, "c" -> WomStringType))
  )

  val nonCoerceableMap = WomMap(WomMapType(WomStringType, WomObjectType), Map(
    WomString("a") -> WomObject(Map.empty),
    WomString("b") -> WomObject(Map.empty),
    WomString("c") -> WomObject(Map.empty))
  )

  "WomObject" should "stringify its value" in {
    abcObject.toWomString shouldEqual "object {a: \"one\", b: \"two\", c: \"three\"}"
  }

  it should "stringify its type" in {
    abcObject.womType.stableName shouldEqual "Object"
  }

  it should "coerce a coerceable map into a WomObject" in {
    WomObjectType.coerceRawValue(coerceableMap) match {
      case Success(v) =>
        v.womType shouldEqual WomObjectType
        v.toWomString shouldEqual abcObject.toWomString
      case Failure(_) => fail("Failed to coerce a map to an object")
    }
  }

  it should "coerce a coerceable typed object into a WomObject" in {
    WomObjectType.coerceRawValue(coerceableTypedObject) match {
      case Success(v) =>
        v.womType shouldEqual WomObjectType
        v.toWomString shouldEqual abcMixedTypedObject.toWomString
      case Failure(_) => fail("Failed to coerce a map to an object")
    }
  }

  it should "NOT successfully coerce a NON coerceable map into a WomObject" in {
    WomObjectType.coerceRawValue(nonCoerceableMap) match {
      case Success(_) => fail("should not have succeeded")
      case Failure(_) => // expected
    }
  }
}
