package wdl.types

import org.scalatest.{FlatSpec, Matchers}
import wdl.WdlExpression
import wdl.expression.NoFunctions
import wom.types._
import wom.values.{WomArray, WomValue}

import scala.util.{Failure, Success}

class WdlArrayTypeSpec extends FlatSpec with Matchers  {
  List(WomStringType, WomArrayType(WomIntegerType), WomPairType(WomIntegerType, WomPairType(WomIntegerType, WomIntegerType)), WomOptionalType(WomStringType)) foreach { desiredMemberType =>
    it should s"be able to construct an empty Array[${desiredMemberType.toDisplayString}] value" in {
      def noLookup(String: String): WomValue = fail("No identifiers should be looked up in this test")

      val desiredArrayType = WomArrayType(desiredMemberType)
      WdlExpression.fromString("[]").evaluate(noLookup, NoFunctions) match {
        case Success(emptyArray @ WomArray(actualArrayType @ WomArrayType(actualMemberType), actualArrayValue)) =>
          actualMemberType should be(WomNothingType)
          actualArrayValue should be(Seq.empty)
          desiredArrayType.isCoerceableFrom(actualArrayType) should be(true)
          desiredArrayType.coerceRawValue(emptyArray) should be(Success(WomArray(desiredArrayType, Seq.empty)))
        case Failure(f) => fail("Unable to create an empty array.", f)
      }
    }

    val arrayType = WomArrayType(desiredMemberType)
    val optionalArrayType = WomArrayType(WomOptionalType(desiredMemberType))
    it should s"be able to coerce $arrayType to $optionalArrayType but not vice versa" in {
      arrayType.isCoerceableFrom(optionalArrayType) should be(false)
      optionalArrayType.isCoerceableFrom(arrayType) should be(true)
    }
  }
}
