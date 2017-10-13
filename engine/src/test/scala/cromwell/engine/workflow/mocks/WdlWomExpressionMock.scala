package cromwell.engine.workflow.mocks

import org.specs2.mock.Mockito
import wdl.{WdlExpression, WdlWomExpression}
import wdl.WdlExpression._
import wdl.expression.WdlFunctions
import wom.values.{WdlString, WdlValue}
import wom.values.{WdlInteger, WdlString, WdlValue}

import scala.util.Success

trait WdlWomExpressionMock extends Mockito {
  val helloStringExpression = {
    val expression = mock[WdlExpression]
    expression.valueString returns "hello"
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WdlValue]]) returns Success(WdlString("hello"))
    expression
  }
  
  def mockStringExpression(value: String): WdlWomExpression = {
    val expression = mock[WdlExpression]
    expression.valueString returns value
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WdlValue]]) returns Success(WdlString(value))
    WdlWomExpression(expression, None)
  }

  def mockIntExpression(value: Int): WdlWomExpression = {
    val expression = mock[WdlExpression]
    expression.valueString returns value.toString
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WdlValue]]) returns Success(WdlInteger(value))
    WdlWomExpression(expression, None)
  }
}
