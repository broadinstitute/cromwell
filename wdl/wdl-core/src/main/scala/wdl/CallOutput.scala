package wdl

import wdl.versioning.WdlVersionSpecifics
import wdl4s.parser.WdlParser.Ast
import wom.types.WomType

final case class CallOutput(call: WdlCall, taskOutput: Output)(implicit implicitVersionSpecifics: WdlVersionSpecifics) extends Output {
  override val wdlVersionSpecifics = implicitVersionSpecifics
  override lazy val requiredExpression: WdlExpression = taskOutput.requiredExpression
  override lazy val ast: Ast = taskOutput.ast
  override lazy val womType: WomType = taskOutput.womType
  override lazy val unqualifiedName: LocallyQualifiedName = taskOutput.unqualifiedName
}
