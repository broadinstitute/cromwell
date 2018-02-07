package cromwell.engine.workflow.mocks

import org.specs2.mock.Mockito
import wdl.{Scope, WdlExpression, WdlWomExpression}
import wdl.WdlExpression._
import wdl.expression.WdlFunctions
import wdl.versioning.NoVersionSpecifics
import wom.values.{WomInteger, WomString, WomValue}

import scala.util.Success

trait WdlWomExpressionMock extends Mockito {
  implicit val wdlVersionSpecifics = NoVersionSpecifics
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
