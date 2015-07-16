package cromwell.binding

import cromwell.binding.values.WdlValue
import org.scalatest.{Matchers, FlatSpec}

class WdlExpressionToFileSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString
  def noLookup(String: String): WdlValue = fail("No identifiers should be looked up in this test")
  class NoFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = fail("No functions should be called in this test")
  }


  "WdlExpressionToFile" should "find no files in an empty expression" in {
    val e: WdlExpression = expr("")
    assert e.anonFiles(noLookup, new NoFunctions()) isEmpty
  }
}
