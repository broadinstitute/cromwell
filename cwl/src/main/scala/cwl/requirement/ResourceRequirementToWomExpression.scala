package cwl.requirement

import cwl.Expression
import shapeless.Poly1
import wom.expression.{ValueAsAnExpression, WomExpression}
import wom.values.{WomInteger, WomString}

object ResourceRequirementToWomExpression extends Poly1 {
  implicit def fromLong: Case.Aux[Long, Set[String] => WomExpression] = at[Long] { l => Function.const(ValueAsAnExpression(WomInteger(l.toInt))) }
  implicit def fromString: Case.Aux[String, Set[String] => WomExpression] = at[String] { s => Function.const(ValueAsAnExpression(WomString(s))) }
  implicit def fromExpression: Case.Aux[Expression, Set[String] => WomExpression] = at[Expression] { e => inputs => 
    cwl.ResourceRequirementExpression(e, inputs)
  }
}
