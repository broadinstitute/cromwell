package wdl4s.wdl

import wdl4s.parser.WdlParser.Ast
import wdl4s.wdl.types.WdlType
import wdl4s.wom.graph.GraphNodePort.GraphNodeOutputPort

object CallOutput {
  def buildOutputPort(callOutput: CallOutput) = {
    GraphNodeOutputPort(callOutput.unqualifiedName, callOutput.wdlType, callOutput.call.womCallNode)
  }
}

case class CallOutput(call: WdlCall, taskOutput: Output) extends Output {
  override lazy val requiredExpression: WdlExpression = taskOutput.requiredExpression
  override lazy val ast: Ast = taskOutput.ast
  override lazy val wdlType: WdlType = taskOutput.wdlType
  override lazy val unqualifiedName: LocallyQualifiedName = taskOutput.unqualifiedName
  
  lazy val toWomOutputPort = CallOutput.buildOutputPort(this)
}
