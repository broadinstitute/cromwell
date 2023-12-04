package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.model.draft3.elements._

object AstToWorkflowBodyElement {
  def astToWorkflowBodyElement(implicit
    astNodeToInputsSectionElement: CheckedAtoB[GenericAstNode, InputsSectionElement],
    astNodeToOutputsSectionElement: CheckedAtoB[GenericAstNode, OutputsSectionElement],
    astNodeToMetaSectionElement: CheckedAtoB[GenericAstNode, MetaSectionElement],
    astNodeToParameterMetaSectionElement: CheckedAtoB[GenericAstNode, ParameterMetaSectionElement],
    astNodeToGraphElement: CheckedAtoB[GenericAstNode, WorkflowGraphElement]
  ): CheckedAtoB[GenericAst, WorkflowBodyElement] = CheckedAtoB.fromCheck { ast: GenericAst =>
    ast.getName match {
      case "Inputs" => astNodeToInputsSectionElement(ast)
      case "Outputs" => astNodeToOutputsSectionElement(ast)
      case "Meta" => astNodeToMetaSectionElement(ast)
      case "ParameterMeta" => astNodeToParameterMetaSectionElement(ast)
      case "Declaration" | "Call" | "Scatter" | "If" => astNodeToGraphElement(ast)
      case other => s"No conversion defined for Ast with name $other to WorkflowBodyElement".invalidNelCheck
    }
  }
}
