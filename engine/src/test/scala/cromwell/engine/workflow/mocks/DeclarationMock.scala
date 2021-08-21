package cromwell.engine.workflow.mocks

import common.mock.MockSugar
import wdl.draft2.model.{Declaration, WdlExpression}
import wom.types.WomType

object DeclarationMock {
  type DeclarationMockType = (String, WomType, WdlExpression)
}

trait DeclarationMock extends MockSugar {
  def mockDeclaration(name: String,
                      womType: WomType,
                      expression: WdlExpression): Declaration = {
    val declaration = mock[Declaration]
    declaration.unqualifiedName returns name
    declaration.expression returns Option(expression)
    declaration.womType returns womType
    declaration
  }
}
