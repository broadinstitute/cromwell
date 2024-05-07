package wdl.transforms.biscayne

import cats.implicits.catsSyntaxValidatedId
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMapTuple2, _}
import wdl.model.draft3.elements.{ExpressionElement, RuntimeAttributesSectionElement}
import wdl.model.draft3.elements.ExpressionElement.{ArrayLiteral, KvPair, PrimitiveLiteralExpressionElement}
import wdl.model.draft3.graph.ExpressionValueConsumer
import wdl.model.draft3.graph.expression.{FileEvaluator, TypeEvaluator, ValueEvaluator}
import wdl.transforms.base.wdlom2wom.TaskDefinitionElementToWomTaskDefinition.{
  eliminateInputDependencies,
  TaskDefinitionElementToWomInputs
}
import wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.transforms.base.wdlom2wom.{
  FileElementToWomBundle,
  FileElementToWomBundleInputs,
  TaskDefinitionElementToWomTaskDefinition,
  WorkflowDefinitionElementToWomWorkflowDefinition
}
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wdl.transforms.biscayne.linking.expression.consumed._
import wdl.transforms.biscayne.linking.expression.files._
import wdl.transforms.biscayne.linking.expression.types._
import wdl.transforms.biscayne.linking.expression.values._
import wom.values.WomInteger
import wom.{RuntimeAttributes, RuntimeAttributesKeys}

package object wdlom2wom {
  val taskDefinitionElementToWomTaskDefinition: CheckedAtoB[TaskDefinitionElementToWomInputs, CallableTaskDefinition] =
    CheckedAtoB.fromErrorOr(convert)
  val workflowDefinitionElementToWomWorkflowDefinition
    : CheckedAtoB[WorkflowDefinitionConvertInputs, WorkflowDefinition] =
    CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomBundle: CheckedAtoB[FileElementToWomBundleInputs, WomBundle] =
    CheckedAtoB.fromCheck(FileElementToWomBundle.convert)

  private def convert(b: TaskDefinitionElementToWomInputs)(implicit
    expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
    fileEvaluator: FileEvaluator[ExpressionElement],
    typeEvaluator: TypeEvaluator[ExpressionElement],
    valueEvaluator: ValueEvaluator[ExpressionElement]
  ): ErrorOr[CallableTaskDefinition] = {
    val a = eliminateInputDependencies(b)(expressionValueConsumer)

    val conversion =
      TaskDefinitionElementToWomTaskDefinition.createTaskGraphAndValidateMetadata(a)(expressionValueConsumer,
                                                                                     fileEvaluator,
                                                                                     typeEvaluator,
                                                                                     valueEvaluator
      ) flatMapN { (taskGraph, _) =>
        val validRuntimeAttributes: ErrorOr[RuntimeAttributes] = a.taskDefinitionElement.runtimeSection match {
          case Some(attributeSection) =>
            TaskDefinitionElementToWomTaskDefinition.createRuntimeAttributes(
              RuntimeAttributesSectionElement(getFinalRuntimeAttributes(attributeSection)),
              taskGraph.linkedGraph
            )(
              expressionValueConsumer,
              fileEvaluator,
              typeEvaluator,
              valueEvaluator
            )
          case None => RuntimeAttributes(Map.empty).validNel
        }

        TaskDefinitionElementToWomTaskDefinition.createCallableTaskDefinition(a, taskGraph, validRuntimeAttributes)(
          expressionValueConsumer,
          fileEvaluator,
          typeEvaluator,
          valueEvaluator
        )
      }

    conversion.contextualizeErrors(s"process task definition '${b.taskDefinitionElement.name}'")
  }

  /**
   * Combine `returnCodes` and `continueOnReturnCode` to be a single attribute. The resulting vector will contain 
   * `continueOnReturnCode` if either `continueOnReturnCode` or `returnCodes` was in the `attributeSection`, it will 
   * never contain `returnCodes`. The value for the `continueOnReturnCode` key in the new vector will be the value 
   * associated with `returnCodes` in the original vector if it exists, else it will be the value associated with 
   * `continueOnReturnCode` in the original vector.
   * @param attributeSection list of all runtime attributes and their values
   * @return A vector of pairs of runtime attribute keys to their respective values
   */
  private def getFinalRuntimeAttributes(attributeSection: RuntimeAttributesSectionElement): Vector[KvPair] = {
    val returnCodesAttribute =
      attributeSection.runtimeAttributes.toList.find(pair => pair.key.equals(RuntimeAttributesKeys.ReturnCodesKey))
    val continueOnReturnCodeAttribute =
      attributeSection.runtimeAttributes.toList.find(pair =>
        pair.key.equals(RuntimeAttributesKeys.ContinueOnReturnCodeKey)
      )

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
        attributeSection.runtimeAttributes.filterNot(attribute =>
          attribute.key.equals(RuntimeAttributesKeys.ContinueOnReturnCodeKey)
        ) ++ Vector(
          KvPair(RuntimeAttributesKeys.ContinueOnReturnCodeKey, returnCodesValue.value)
        )
      case _ => attributeSection.runtimeAttributes
    }

    finalAttributes.filterNot(attribute => attribute.key.equals(RuntimeAttributesKeys.ReturnCodesKey))
  }
}
