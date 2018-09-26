package wdl.transforms.biscayne

import common.transforms.CheckedAtoB
import wdl.transforms.base.wdlom2wom.TaskDefinitionElementToWomTaskDefinition.TaskDefinitionElementToWomInputs
import wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.transforms.base.wdlom2wom.{FileElementToWomBundle, FileElementToWomBundleInputs, TaskDefinitionElementToWomTaskDefinition, WorkflowDefinitionElementToWomWorkflowDefinition}
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wdl.transforms.biscayne.linking.expression.consumed._
import wdl.transforms.biscayne.linking.expression.files._
import wdl.transforms.biscayne.linking.expression.types._
import wdl.transforms.biscayne.linking.expression.values._

package object wdlom2wom {
  val taskDefinitionElementToWomTaskDefinition: CheckedAtoB[TaskDefinitionElementToWomInputs, CallableTaskDefinition] = CheckedAtoB.fromErrorOr(TaskDefinitionElementToWomTaskDefinition.convert)
  val workflowDefinitionElementToWomWorkflowDefinition: CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] = CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomBundle: CheckedAtoB[FileElementToWomBundleInputs, WomBundle] = CheckedAtoB.fromCheck(FileElementToWomBundle.convert)
}
