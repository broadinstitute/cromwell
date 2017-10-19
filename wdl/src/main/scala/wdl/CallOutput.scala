package wdl

import wdl4s.parser.WdlParser.Ast
import wom.types.WdlType

final case class CallOutput(call: WdlCall, taskOutput: Output) extends Output {
  override lazy val requiredExpression: WdlExpression = taskOutput.requiredExpression
  override lazy val ast: Ast = taskOutput.ast
  override lazy val wdlType: WdlType = taskOutput.wdlType
  override lazy val unqualifiedName: LocallyQualifiedName = taskOutput.unqualifiedName
}
