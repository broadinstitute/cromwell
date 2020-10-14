package wdl.values

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft2.model.expression.NoFunctions
import wdl.draft2.model.{NoLookup, WdlExpression}
import wom.types.{WomArrayType, WomIntegerType}
import wom.values.WomArray

import scala.util.Success

class LargeArraySpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  it should s"not have a stack overflow when parsing an array of 100,000 elements" in {
    val expr = WdlExpression.fromString(s"[${(1 to 100000).mkString(",")}]")
    expr.evaluate(NoLookup, NoFunctions) match {
      case Success(a: WomArray) =>
        a.value.size shouldEqual 100000
        a.womType shouldEqual WomArrayType(WomIntegerType)
      case _ => fail("Expecting to parse as an Array[Int]")
    }
  }
}
