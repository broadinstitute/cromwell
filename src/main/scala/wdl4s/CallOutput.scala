package wdl4s

import wdl4s.parser.WdlParser.Ast
import wdl4s.types.WdlType

case class CallOutput(call: Call, taskOutput: Output) extends Output {
  override lazy val requiredExpression: WdlExpression = taskOutput.requiredExpression
  override lazy val ast: Ast = taskOutput.ast
  override lazy val wdlType: WdlType = taskOutput.wdlType
  override lazy val unqualifiedName: LocallyQualifiedName = taskOutput.unqualifiedName
}
