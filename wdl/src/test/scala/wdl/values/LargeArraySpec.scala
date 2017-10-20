package wdl.values

import org.scalatest.{FlatSpec, Matchers}
import wdl.expression.NoFunctions
import wdl.{NoLookup, WdlExpression}
import wom.types.{WomArrayType, WomIntegerType}
import wom.values.WomArray

import scala.util.Success

class LargeArraySpec extends FlatSpec with Matchers {
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
