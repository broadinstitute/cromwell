package wdl.draft3.transforms

import common.transforms.CheckedAtoB
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wom.callable.WorkflowDefinition
import wom.executable.WomBundle

package object wdlom2wom {
  val workflowDefinitionElementToWomWorkflowDefinition: CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] = CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomBundle: CheckedAtoB[FileElementAndImportResolvers, WomBundle] = CheckedAtoB.fromCheck(FileElementToWomBundle.convert)
}
