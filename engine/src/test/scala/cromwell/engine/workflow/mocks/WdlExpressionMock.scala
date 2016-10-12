package cromwell.engine.workflow.mocks

import org.specs2.mock.Mockito
import wdl4s.WdlExpression
import wdl4s.WdlExpression._
import wdl4s.expression.WdlFunctions
import wdl4s.values.{WdlInteger, WdlString, WdlValue}

import scala.util.Success

trait WdlExpressionMock extends Mockito {
  val helloStringExpression = {
    val expression = mock[WdlExpression]
    expression.valueString returns "hello"
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WdlValue]]) returns Success(WdlString("hello"))
    expression
  }
  
  def mockStringExpression(value: String) = {
    val expression = mock[WdlExpression]
    expression.valueString returns value
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WdlValue]]) returns Success(WdlString(value))
    expression
  }

  def mockIntExpression(value: Int) = {
    val expression = mock[WdlExpression]
    expression.valueString returns value.toString
    expression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WdlValue]]) returns Success(WdlInteger(value))
    expression
  }
}
