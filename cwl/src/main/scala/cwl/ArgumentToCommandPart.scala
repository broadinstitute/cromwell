package cwl

import cwl.command.StringCommandPart
import shapeless._
import wom.CommandPart

object ArgumentToCommandPart extends Poly1 {
  type ExpressionToCommandPart = ExpressionLib => CommandPart
  implicit def caseStringOrExpression: Case.Aux[StringOrExpression, ExpressionToCommandPart] = at {
    _.fold(this)
  }

  implicit def caseExpression: Case.Aux[Expression, ExpressionToCommandPart] = at {
    CwlExpressionCommandPart.apply
  }

  implicit def caseString: Case.Aux[String, ExpressionToCommandPart] = at {
    string =>
      _ =>
        StringCommandPart(string)
  }

  implicit def caseCommandLineBinding: Case.Aux[ArgumentCommandLineBinding, ExpressionToCommandPart] = at {
    ArgumentCommandLineBindingCommandPart.apply
  }
}
