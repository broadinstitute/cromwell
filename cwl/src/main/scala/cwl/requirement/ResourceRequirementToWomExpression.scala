package cwl.requirement

import cwl.{Expression, ExpressionLib}
import shapeless.Poly1
import wom.expression.{ValueAsAnExpression, WomExpression}
import wom.values.{WomLong, WomString}

object ResourceRequirementToWomExpression extends Poly1 {
  type ResourceRequirementStringSetToWomExpression = (Set[String], ExpressionLib) => WomExpression
  implicit def fromLong: Case.Aux[Long, ResourceRequirementStringSetToWomExpression] = at[Long] { l => (_, _) => ValueAsAnExpression(WomLong(l)) }
  implicit def fromString: Case.Aux[String, ResourceRequirementStringSetToWomExpression] = at[String] { s =>  (_, _) => ValueAsAnExpression(WomString(s)) }
  implicit def fromExpression: Case.Aux[Expression, ResourceRequirementStringSetToWomExpression] = at[Expression] { e => (inputs, expressionLib) =>
    cwl.ECMAScriptWomExpression(e, inputs, expressionLib)
  }
}
