package wdl

import wdl4s.parser.WdlParser.Ast
import wom.types.WomType

final case class CallOutput(call: WdlCall, taskOutput: Output) extends Output {
  override lazy val requiredExpression: WdlExpression = taskOutput.requiredExpression
  override lazy val ast: Ast = taskOutput.ast
  override lazy val womType: WomType = taskOutput.womType
  override lazy val unqualifiedName: LocallyQualifiedName = taskOutput.unqualifiedName
}
