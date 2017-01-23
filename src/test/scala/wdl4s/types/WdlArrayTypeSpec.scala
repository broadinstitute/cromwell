package wdl4s.types

import wdl4s.values.{WdlArray, WdlInteger, WdlOptionalValue, WdlString, WdlValue}
import wdl4s.parser.WdlParser.SyntaxError
import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsArray, JsNumber}
import wdl4s.WdlExpression
import wdl4s.expression.NoFunctions

import scala.util.{Failure, Success}

class WdlArrayTypeSpec extends FlatSpec with Matchers  {
  val intArray = WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(1), WdlInteger(2), WdlInteger(3)))
  "WdlArray" should "stringify its value" in {
    intArray.toWdlString shouldEqual "[1, 2, 3]"
  }
  it should "tsv serialize an Array[Array[String]]" in {
    val nestedArray = WdlArray(
      WdlArrayType(WdlArrayType(WdlStringType)),
      Seq(
        WdlArray(WdlArrayType(WdlStringType), Seq("a", "b").map(WdlString)),
        WdlArray(WdlArrayType(WdlStringType), Seq("c", "d").map(WdlString))
      )
    )
    nestedArray.tsvSerialize shouldEqual Success("a\tb\nc\td")
  }
  it should "fail to TSV serialize an Array[Array[Array[String]]]" in {
    try {
      val array = WdlArray(WdlArrayType(WdlStringType), Seq("a", "b").map(WdlString))
      val nestedArray = WdlArray(
        WdlArrayType(WdlArrayType(WdlStringType)),
        Seq(
          WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), Seq(array, array)),
          WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), Seq(array, array))
        )
      )
      nestedArray.tsvSerialize
      fail("should not have succeeded")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  "WdlArrayType" should "coerce Seq(1,2,3) into a WdlArray" in {
    WdlArrayType(WdlIntegerType).coerceRawValue(Seq(1,2,3)) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing array: $f")
    }
  }
  it should "coerce a JsArray into a WdlArray" in {
    WdlArrayType(WdlIntegerType).coerceRawValue(JsArray(JsNumber(1), JsNumber(2), JsNumber(3))) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing JsArray: $f")
    }
  }
  it should "coerce single values into one-element arrays" in {
    WdlArrayType(WdlStringType).coerceRawValue(WdlString("edamame is tasty")) match {
      case Success(array) => array shouldEqual WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("edamame is tasty")))
      case Failure(f) => fail("exception coercing single value to array", f)
    }
  }
  it should "stringify its type" in {
    intArray.wdlType.toWdlString shouldEqual "Array[Int]"
  }
  it should "convert WDL source code to WdlArray" in {
    WdlArrayType(WdlIntegerType).fromWdlString("[1,2,3]") shouldEqual intArray
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST" in {
    try {
      WdlArrayType(WdlIntegerType).fromWdlString("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (2)" in {
    try {
      WdlArrayType(WdlIntegerType).fromWdlString("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (3)" in {
    try {
      WdlArrayType(WdlIntegerType).fromWdlString("[1,stdout()]")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (4)" in {
    try {
      WdlArrayType(WdlIntegerType).fromWdlString("[1,var]")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "detect invalid array construction if there are uncoerceable types" in {
    try {
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("foo"), WdlOptionalValue(WdlStringType, None)))
      fail("Invalid array initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "detect invalid array construction if type does not match the input array type" in {
    try {
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlInteger(2)))
      fail("Invalid array initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }

  List(WdlStringType, WdlArrayType(WdlIntegerType), WdlPairType(WdlIntegerType, WdlPairType(WdlIntegerType, WdlIntegerType)), WdlOptionalType(WdlStringType)) foreach { desiredMemberType =>
    it should s"be able to construct an empty Array[${desiredMemberType.toWdlString}] value" in {
      def noLookup(String: String): WdlValue = fail("No identifiers should be looked up in this test")

      val desiredArrayType = WdlArrayType(desiredMemberType)
      WdlExpression.fromString("[]").evaluate(noLookup, NoFunctions) match {
        case Success(emptyArray @ WdlArray(actualArrayType @ WdlArrayType(actualMemberType), actualArrayValue)) =>
          actualMemberType should be(WdlAnyType)
          actualArrayValue should be(Seq.empty)
          desiredArrayType.isCoerceableFrom(actualArrayType) should be(true)
          desiredArrayType.coerceRawValue(emptyArray) should be(Success(WdlArray(desiredArrayType, Seq.empty)))
        case Failure(f) => fail("Unable to create an empty array.", f)
      }

    }
  }
}
