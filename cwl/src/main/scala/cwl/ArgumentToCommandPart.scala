package cwl

import cwl.command.StringCommandPart
import shapeless._
import wom.CommandPart

object ArgumentToCommandPart extends Poly1 {
  implicit def caseStringOrExpression: Case.Aux[StringOrExpression, CommandPart] = at {
    _.fold(this)
  }

  implicit def caseExpression: Case.Aux[Expression, CommandPart] = at {
    CwlExpressionCommandPart.apply
  }

  implicit def caseString: Case.Aux[String, CommandPart] = at {
    StringCommandPart.apply
  }

  implicit def caseCommandLineBinding: Case.Aux[ArgumentCommandLineBinding, CommandPart] = at {
    ArgumentCommandLineBindingCommandPart.apply
  }
}
