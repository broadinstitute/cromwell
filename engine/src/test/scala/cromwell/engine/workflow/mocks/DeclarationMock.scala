package cromwell.engine.workflow.mocks

import org.specs2.mock.Mockito
import wdl.draft2.model.{Declaration, WdlExpression}
import wom.types.WomType

object DeclarationMock {
  type DeclarationMockType = (String, WomType, WdlExpression)
}

trait DeclarationMock extends Mockito {
  def mockDeclaration(name: String,
                      womType: WomType,
                      expression: WdlExpression) = {
    val declaration = mock[Declaration]
    declaration.unqualifiedName returns name
    declaration.expression returns Option(expression)
    declaration.womType returns womType
    declaration
  }
}
