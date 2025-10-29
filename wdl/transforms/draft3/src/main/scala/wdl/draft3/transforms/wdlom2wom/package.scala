package wdl.draft3.transforms

import common.transforms.CheckedAtoB
import wdl.transforms.base.wdlom2wom.TaskDefinitionElementToWomTaskDefinition.TaskDefinitionElementToWomInputs
import wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.transforms.base.wdlom2wom.{
  FileElementToWomBundle,
  FileElementToWomBundleInputs,
  TaskDefinitionElementToWomTaskDefinition,
  WorkflowDefinitionElementToWomWorkflowDefinition
}
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wdl.draft3.transforms.linking.expression.consumed._
import wdl.draft3.transforms.linking.expression.files._
import wdl.draft3.transforms.linking.expression.types._
import wdl.draft3.transforms.linking.expression.values._
import wdl.model.draft3.elements.RuntimeAttributesSectionElement
import wom.RuntimeAttributesKeys

package object wdlom2wom {
  val taskDefinitionElementToWomTaskDefinition: CheckedAtoB[TaskDefinitionElementToWomInputs, CallableTaskDefinition] =
    CheckedAtoB.fromErrorOr(a => TaskDefinitionElementToWomTaskDefinition.convert(a, removeFutureAttrs))
  val workflowDefinitionElementToWomWorkflowDefinition
    : CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] =
    CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomBundle: CheckedAtoB[FileElementToWomBundleInputs, WomBundle] =
    CheckedAtoB.fromCheck(FileElementToWomBundle.convert)

  /**
   * Remove attributes that are not supported in WDL 1.0
   *   - `container` supercedes `docker` starting with WDL 1.1.
   *   - `gpu` controls a GPU requirement check starting in WDL 1.1
   */
  private def removeFutureAttrs(
    attributeSection: Option[RuntimeAttributesSectionElement]
  ): Option[RuntimeAttributesSectionElement] =
    attributeSection.map(_.runtimeAttributes) map { originalAttrs =>
      val finalAttributes = originalAttrs.filterNot(pair =>
        List(RuntimeAttributesKeys.ContainerKey, RuntimeAttributesKeys.GpuRequiredKey).contains(pair.key)
      )
      RuntimeAttributesSectionElement(finalAttributes)
    }
}
