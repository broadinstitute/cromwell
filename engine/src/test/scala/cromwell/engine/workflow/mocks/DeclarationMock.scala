package cromwell.engine.workflow.mocks

import org.specs2.mock.Mockito
import wdl.{Declaration, WdlExpression}
import wom.types.WdlType

object DeclarationMock {
  type DeclarationMockType = (String, WdlType, WdlExpression)
}

trait DeclarationMock extends Mockito {
  def mockDeclaration(name: String,
                      wdlType: WdlType,
                      expression: WdlExpression) = {
    val declaration = mock[Declaration]
    declaration.unqualifiedName returns name
    declaration.expression returns Option(expression)
    declaration.wdlType returns wdlType
    declaration
  }
}
