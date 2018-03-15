package wdl.draft3.transforms

import common.transforms.CheckedAtoB
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.model.draft3.elements.{CommandPartElement, TaskDefinitionElement}
import wom.{CommandPart}
import wom.callable.{TaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle

package object wdlom2wom {
  val commandPartElementToWomCommandPart: CheckedAtoB[CommandPartElement, CommandPart] = CheckedAtoB.fromErrorOr(CommandPartElementToWomCommandPart.convert)
  val taskDefinitionElementToWomTaskDefinition: CheckedAtoB[TaskDefinitionElement, TaskDefinition] = CheckedAtoB.fromErrorOr(TaskDefinitionElementToWomTaskDefinition.convert)
  val workflowDefinitionElementToWomWorkflowDefinition: CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] = CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomBundle: CheckedAtoB[FileElementToWomBundleInputs, WomBundle] = CheckedAtoB.fromCheck(FileElementToWomBundle.convert)
}
