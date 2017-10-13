package wom.types

import org.scalatest.{FlatSpec, Matchers}
import wom.values.{WdlMap, WdlObject, WdlString}

import scala.util.{Failure, Success}

class WdlObjectTypeSpec extends FlatSpec with Matchers {
  val abcObject = WdlObject(Map(
    "a" -> WdlString("one"),
    "b" -> WdlString("two"),
    "c" -> WdlString("three")
  ))

  val coerceableMap = WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
    WdlString("a") -> WdlString("one"),
    WdlString("b") -> WdlString("two"),
    WdlString("c") -> WdlString("three"))
  )

  val nonCoerceableMap = WdlMap(WdlMapType(WdlStringType, WdlObjectType), Map(
    WdlString("a") -> WdlObject(Map.empty),
    WdlString("b") -> WdlObject(Map.empty),
    WdlString("c") -> WdlObject(Map.empty))
  )
  
  "WdlObject" should "stringify its value" in {
    abcObject.toWdlString shouldEqual "object {a: \"one\", b: \"two\", c: \"three\"}"
  }

  it should "stringify its type" in {
    abcObject.wdlType.toWdlString shouldEqual "Object"
  }

  it should "coerce a coerceable map into a WdlObject" in {
    WdlObjectType.coerceRawValue(coerceableMap) match {
      case Success(v) => 
        v.wdlType shouldEqual WdlObjectType
        v.toWdlString shouldEqual abcObject.toWdlString
      case Failure(_) => fail("Failed to coerce a map to an object")
    }
  }

  it should "NOT successfully coerce a NON coerceable map into a WdlObject" in {
    WdlObjectType.coerceRawValue(nonCoerceableMap) match {
      case Success(_) => fail("should not have succeeded")
      case Failure(_) => // expected
    }
  }
}
