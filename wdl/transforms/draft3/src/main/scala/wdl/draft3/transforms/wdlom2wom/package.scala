package wdl.draft3.transforms

import common.transforms.CheckedAtoB

package object wdlom2wom {
  implicit val checkedWorkflowDefinitionElementToWomWorkflowDefinition = CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  implicit val fileElementToWomExecutable = CheckedAtoB.fromErrorOr(FileElementToWomExecutable.convert)
}
