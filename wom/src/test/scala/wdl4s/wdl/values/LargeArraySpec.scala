package wdl4s.wdl.values

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.expression.NoFunctions
import wdl4s.wdl.types.{WdlArrayType, WdlIntegerType}
import wdl4s.wdl.{NoLookup, WdlExpression}
import scala.util.Success

class LargeArraySpec extends FlatSpec with Matchers {
  it should s"not have a stack overflow when parsing an array of 100,000 elements" in {
    val expr = WdlExpression.fromString(s"[${(1 to 100000).mkString(",")}]")
    expr.evaluate(NoLookup, NoFunctions) match {
      case Success(a: WdlArray) =>
        a.value.size shouldEqual 100000
        a.wdlType shouldEqual WdlArrayType(WdlIntegerType)
      case _ => fail("Expecting to parse as an Array[Int]")
    }
  }
}
