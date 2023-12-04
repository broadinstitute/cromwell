package wdl.transforms.base.ast2wdlom

import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.model.draft3.elements.{FileBodyElement, StructElement, TaskDefinitionElement, WorkflowDefinitionElement}

object AstToFileBodyElement {

  def astToFileBodyElement(implicit
    workflowConverter: CheckedAtoB[GenericAst, WorkflowDefinitionElement],
    taskConverter: CheckedAtoB[GenericAst, TaskDefinitionElement],
    structConverter: CheckedAtoB[GenericAst, StructElement]
  ): CheckedAtoB[GenericAst, FileBodyElement] =
    CheckedAtoB.fromCheck(convert(workflowConverter, taskConverter, structConverter))

  def convert(workflowConverter: CheckedAtoB[GenericAst, WorkflowDefinitionElement],
              taskConverter: CheckedAtoB[GenericAst, TaskDefinitionElement],
              structConverter: CheckedAtoB[GenericAst, StructElement]
  )(ast: GenericAst): Checked[FileBodyElement] = ast.getName match {
    case "Workflow" => workflowConverter.run(ast)
    case "Task" => taskConverter.run(ast)
    case "Struct" => structConverter.run(ast)
    case other => s"No conversion defined for Ast with name $other".invalidNelCheck
  }
}
