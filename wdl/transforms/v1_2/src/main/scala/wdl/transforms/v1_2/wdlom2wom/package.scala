package wdl.transforms.v1_2

import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.ExpressionElement.{ArrayLiteral, KvPair, PrimitiveLiteralExpressionElement}
import wdl.model.draft3.elements.RuntimeAttributesSectionElement
import wdl.transforms.base.wdlom2wom.TaskDefinitionElementToWomTaskDefinition.TaskDefinitionElementToWomInputs
import wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.transforms.base.wdlom2wom.{
  FileElementToWomBundle,
  FileElementToWomBundleInputs,
  TaskDefinitionElementToWomTaskDefinition,
  WorkflowDefinitionElementToWomWorkflowDefinition
}
import wdl.transforms.v1_2.linking.expression.consumed._
import wdl.transforms.v1_2.linking.expression.files._
import wdl.transforms.v1_2.linking.expression.types._
import wdl.transforms.v1_2.linking.expression.values._
import wom.RuntimeAttributesKeys
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wom.values.WomInteger

package object wdlom2wom {
  // TODO: WDL 1.2 overrides requirements and hints in the inputs JSON with `wf.task.requirements.<key>` and
  //  `wf.task.hints.<key>`, and they always supersede the values in the WDL. Cromwell's runtime override input is
  //  named `runtime` (`RuntimeOverrideInputDefinition`), so 1.2 users override them with `wf.task.runtime`.
  //  Suggestion: make the override input name(s) configurable per language version in `RuntimeOverrideInputDefinition`
  //  and accept both `requirements` and `hints` for WDL 1.2.
  // TODO: Per the spec, requirements have spec-defined defaults (e.g. 1 cpu, 2 GiB memory) and engines "should not
  //  provide additional mechanisms to set default values". Cromwell still applies `default_runtime_attributes` from
  //  workflow options and backend config. Suggestion: decide whether WDL 1.2 tasks should ignore those defaults.
  // TODO: Call caching hashes the runtime attributes that backends mark as cache-worthy, so hint-sourced attributes
  //  (e.g. `preemptible`) are treated like their WDL 1.1 runtime equivalents. Suggestion: confirm that hints that do
  //  not change outputs (e.g. `short_task`, `localization_optional`) never affect call caching.
  // TODO: WDL 1.2 changes two coercion rules: a Map or Object may be coerced to a Struct even when it has extra keys
  //  (they are ignored), and a Struct may be coerced to another Struct with compatible members. Coercion is implemented
  //  in `wom.types.WomCompositeType`, which is shared by all WDL versions, and ignoring extra keys would silently change
  //  WDL 1.0/1.1 behaviour. Suggestion: add a version-specific coercion mode to `WomCompositeType` (or coerce in the
  //  WDL 1.2 value evaluators) and decide whether to warn when extra keys are dropped.
  val taskDefinitionElementToWomTaskDefinition: CheckedAtoB[TaskDefinitionElementToWomInputs, CallableTaskDefinition] =
    CheckedAtoB.fromErrorOr(a =>
      TaskDefinitionElementToWomTaskDefinition.convert(a,
                                                       createRuntimeOverrideInputs = true,
                                                       combineReturnCodeRuntimeAttributes
      )
    )
  val workflowDefinitionElementToWomWorkflowDefinition
    : CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] =
    CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomBundle: CheckedAtoB[FileElementToWomBundleInputs, WomBundle] =
    CheckedAtoB.fromCheck(FileElementToWomBundle.convert)

  /**
   * Combine `returnCodes` and `continueOnReturnCode` to be a single attribute. The resulting vector will contain 
   * `continueOnReturnCode` if either `continueOnReturnCode` or `returnCodes` was in the `attributeSection`, it will 
   * never contain `returnCodes`. The value for the `continueOnReturnCode` key in the new vector will be the value 
   * associated with `returnCodes` in the original vector if it exists, else it will be the value associated with 
   * `continueOnReturnCode` in the original vector.
   * @param attributeSection list of all runtime attributes and their values
   * @return A vector of pairs of runtime attribute keys to their respective values
   */
  private def combineReturnCodeRuntimeAttributes(
    attributeSection: Option[RuntimeAttributesSectionElement]
  ): Option[RuntimeAttributesSectionElement] =
    attributeSection.map(_.runtimeAttributes) map { originalAttrs =>
      val returnCodesAttribute =
        originalAttrs.toList.find(pair => pair.key.equals(RuntimeAttributesKeys.ReturnCodesKey))
      val continueOnReturnCodeAttribute =
        originalAttrs.toList.find(pair => pair.key.equals(RuntimeAttributesKeys.ContinueOnReturnCodeKey))

      val returnCodesNotUnique = (returnCodesAttribute, continueOnReturnCodeAttribute) match {
        case (Some(returnCodesValue), Some(continueOnReturnCodeValue)) =>
          returnCodesValue.value
            .equals(continueOnReturnCodeValue.value) || returnCodesValue.value.equals(
            ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(0))))
          )
        case _ => false
      }

      val finalAttributes = (returnCodesAttribute, returnCodesNotUnique) match {
        case (Some(returnCodesValue), false) =>
          originalAttrs.filterNot(attribute =>
            attribute.key.equals(RuntimeAttributesKeys.ContinueOnReturnCodeKey)
          ) ++ Vector(
            KvPair(RuntimeAttributesKeys.ContinueOnReturnCodeKey, returnCodesValue.value)
          )
        case _ => originalAttrs
      }

      RuntimeAttributesSectionElement(
        finalAttributes.filterNot(attribute => attribute.key.equals(RuntimeAttributesKeys.ReturnCodesKey))
      )
    }
}
