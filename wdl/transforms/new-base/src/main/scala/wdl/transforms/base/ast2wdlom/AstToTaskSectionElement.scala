package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.model.draft3.elements._

object AstToTaskSectionElement {
  def astToTaskSectionElement(implicit
    astNodeToInputsSectionElement: CheckedAtoB[GenericAstNode, InputsSectionElement],
    astNodeToOutputsSectionElement: CheckedAtoB[GenericAstNode, OutputsSectionElement],
    astNodeToDeclarationContent: CheckedAtoB[GenericAstNode, DeclarationContent],
    astNodeToRuntimeAttributesSectionElement: CheckedAtoB[GenericAstNode, RuntimeAttributesSectionElement],
    astNodeToCommandSectionElement: CheckedAtoB[GenericAstNode, CommandSectionElement],
    astNodeToMetaSectionElement: CheckedAtoB[GenericAstNode, MetaSectionElement],
    astNodeToParameterMetaSectionElement: CheckedAtoB[GenericAstNode, ParameterMetaSectionElement]
  ): CheckedAtoB[GenericAst, TaskSectionElement] = CheckedAtoB.fromCheck { a: GenericAst =>
    a.getName match {
      case "Inputs" => astNodeToInputsSectionElement.run(a)
      case "Outputs" => astNodeToOutputsSectionElement(a)
      case "Declaration" => astNodeToDeclarationContent(a).map(IntermediateValueDeclarationElement.fromContent)
      case "Runtime" => astNodeToRuntimeAttributesSectionElement(a)
      case "RawCommand" => astNodeToCommandSectionElement(a)
      case "Meta" => astNodeToMetaSectionElement(a)
      case "ParameterMeta" => astNodeToParameterMetaSectionElement(a)
      case other => s"No conversion defined for Ast with name $other to TaskBodyElement".invalidNelCheck
    }
  }
}
