package cwl

import cwl.command.StringCommandPart
import shapeless._
import wom.CommandPart

object ArgumentToCommandPart extends Poly1 {
  type MakeCommandPart = (Boolean, ExpressionLib) => CommandPart
  implicit def caseStringOrExpression: Case.Aux[StringOrExpression, MakeCommandPart] = at {
    _.fold(this)
  }

  implicit def caseExpression: Case.Aux[Expression, MakeCommandPart] = at {
    CwlExpressionCommandPart.apply
  }

  implicit def caseString: Case.Aux[String, MakeCommandPart] = at {
    string =>
      (_, _) =>
        StringCommandPart(string)
  }

  implicit def caseCommandLineBinding: Case.Aux[ArgumentCommandLineBinding, MakeCommandPart] = at {
    ArgumentCommandLineBindingCommandPart.apply
  }
}
