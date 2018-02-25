package cromwell.engine.workflow.mocks

import org.specs2.mock.Mockito
import wdl.draft2.model.{Scope, WdlExpression, WdlWomExpression}
import wdl.draft2.model.WdlExpression._
import wdl.draft2.model.expression.WdlFunctions
import wom.values.{WomInteger, WomString, WomValue}

import scala.util.Success

trait WdlWomExpressionMock extends Mockito {
  val helloStringExpression = {
    val expression = mock[WdlExpression]
    expression.valueString returns "hello"
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WomValue]]) returns Success(WomString("hello"))
    expression
  }
  
  def mockStringExpression(value: String): WdlWomExpression = {
    val expression = mock[WdlExpression]
    expression.valueString returns value
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WomValue]]) returns Success(WomString(value))

    val mockScope = mock[Scope]
    WdlWomExpression(expression, mockScope)
  }

  def mockIntExpression(value: Int): WdlWomExpression = {
    val expression = mock[WdlExpression]
    expression.valueString returns value.toString
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WomValue]]) returns Success(WomInteger(value))

    val mockScope = mock[Scope]
    WdlWomExpression(expression, mockScope)
  }
}
