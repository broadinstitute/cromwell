package wdl.transforms.cascades

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
import wdl.transforms.cascades.linking.expression.consumed._
import wdl.transforms.cascades.linking.expression.files._
import wdl.transforms.cascades.linking.expression.types._
import wdl.transforms.cascades.linking.expression.values._
import wom.RuntimeAttributesKeys
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wom.values.WomInteger

package object wdlom2wom {
  val taskDefinitionElementToWomTaskDefinition: CheckedAtoB[TaskDefinitionElementToWomInputs, CallableTaskDefinition] =
    CheckedAtoB.fromErrorOr(a =>
      TaskDefinitionElementToWomTaskDefinition.convert(a, combineReturnCodeRuntimeAttributes)
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
