package wdl.types

import org.scalatest.{FlatSpec, Matchers}
import wdl.WdlExpression
import wdl.expression.NoFunctions
import wom.types._
import wom.values.{WdlArray, WdlValue}

import scala.util.{Failure, Success}

class WdlArrayTypeSpec extends FlatSpec with Matchers  {
  List(WdlStringType, WdlArrayType(WdlIntegerType), WdlPairType(WdlIntegerType, WdlPairType(WdlIntegerType, WdlIntegerType)), WdlOptionalType(WdlStringType)) foreach { desiredMemberType =>
    it should s"be able to construct an empty Array[${desiredMemberType.toWdlString}] value" in {
      def noLookup(String: String): WdlValue = fail("No identifiers should be looked up in this test")

      val desiredArrayType = WdlArrayType(desiredMemberType)
      WdlExpression.fromString("[]").evaluate(noLookup, NoFunctions) match {
        case Success(emptyArray @ WdlArray(actualArrayType @ WdlArrayType(actualMemberType), actualArrayValue)) =>
          actualMemberType should be(WdlNothingType)
          actualArrayValue should be(Seq.empty)
          desiredArrayType.isCoerceableFrom(actualArrayType) should be(true)
          desiredArrayType.coerceRawValue(emptyArray) should be(Success(WdlArray(desiredArrayType, Seq.empty)))
        case Failure(f) => fail("Unable to create an empty array.", f)
      }
    }

    val arrayType = WdlArrayType(desiredMemberType)
    val optionalArrayType = WdlArrayType(WdlOptionalType(desiredMemberType))
    it should s"be able to coerce $arrayType to $optionalArrayType but not vice versa" in {
      arrayType.isCoerceableFrom(optionalArrayType) should be(false)
      optionalArrayType.isCoerceableFrom(arrayType) should be(true)
    }
  }
}
