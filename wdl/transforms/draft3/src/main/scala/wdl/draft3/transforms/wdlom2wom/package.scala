package wdl.draft3.transforms

import common.transforms.CheckedAtoB
import wom.executable.WomBundle

package object wdlom2wom {
  val workflowDefinitionElementToWomWorkflowDefinition = CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomBundle: CheckedAtoB[FileElementAndImportResolvers, WomBundle] = CheckedAtoB.fromCheck(FileElementToWomBundle.convert)
}
