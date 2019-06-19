package wdl.types

import org.scalatest.{FlatSpec, Matchers}
import wdl.draft2.model.WdlExpression
import wdl.draft2.model.expression.NoFunctions
import wom.types._
import wom.values.{WomArray, WomValue}

import scala.util.{Failure, Success}

class WomArrayTypeSpec extends FlatSpec with Matchers  {

  behavior of "WomArrayType"

  List(WomStringType, WomArrayType(WomIntegerType), WomPairType(WomIntegerType, WomPairType(WomIntegerType, WomIntegerType)), WomOptionalType(WomStringType)) foreach { desiredMemberType =>
    it should s"be able to construct an empty Array[${desiredMemberType.stableName}] value" in {

      val desiredArrayType = WomArrayType(desiredMemberType)
      WdlExpression.fromString("[]").evaluate(noLookup, NoFunctions) match {
        case Success(emptyArray @ WomArray(actualArrayType @ WomMaybeEmptyArrayType(actualMemberType), actualArrayValue)) =>
          actualMemberType should be(WomNothingType)
          actualArrayValue should be(Seq.empty)
          desiredArrayType.isCoerceableFrom(actualArrayType) should be(true)
          desiredArrayType.coerceRawValue(emptyArray) should be(Success(WomArray(desiredArrayType, Seq.empty)))
        case Success(WomArray(WomNonEmptyArrayType(_), _)) => fail("Empty arrays should not be created with an Array[_]+ type")
        case Success(other) => fail(s"Array literal somehow got evaluated as a ${other.womType} type?!?")
        case Failure(f) => fail("Unable to create an empty array.", f)
      }
    }

    val arrayType = WomArrayType(desiredMemberType)
    val optionalArrayType = WomArrayType(WomOptionalType(desiredMemberType))
    it should s"be able to coerce ${arrayType.stableName} to ${optionalArrayType.stableName} and vice versa" in {

      // If the inner type is already an optional, we can use the flatten coercion to go backwards (eg from Array[String??] to Array[String?])
      val backwardsCoerceable = arrayType.memberType.isInstanceOf[WomOptionalType]
      arrayType.isCoerceableFrom(optionalArrayType) should be(backwardsCoerceable)
      optionalArrayType.isCoerceableFrom(arrayType) should be(true)
    }
  }

  val nonEmptyLiteral = "[1, 2, 3]"

  List(WomNonEmptyArrayType(WomIntegerType), WomOptionalType(WomNonEmptyArrayType(WomIntegerType))) foreach { targetType =>
    it should s"be able to coerce an array literal into ${targetType.stableName}" in {

      val evaluatedLiteral: WomValue = WdlExpression.fromString(nonEmptyLiteral).evaluate(noLookup, NoFunctions).getOrElse(fail(s"Unable to evaluate non-empty literal $nonEmptyLiteral"))

      targetType.coerceRawValue(evaluatedLiteral) match {
        case Success(womValue) => womValue.womType should be(targetType)
        case Failure(e) => fail(s"Unable to coerce $evaluatedLiteral (${evaluatedLiteral.womType.stableName}) into ${targetType.stableName}", e)
      }
    }
  }

  def noLookup(String: String): WomValue = fail("No identifiers should be looked up in this test")
}
