package wdl.transforms.cascades

import cats.implicits.{catsSyntaxTuple2Semigroupal, catsSyntaxValidatedId, toTraverseOps}
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMapTuple2, _}
import wdl.model.draft3.elements.ExpressionElement.{ArrayLiteral, KvPair, PrimitiveLiteralExpressionElement}
import wdl.model.draft3.elements.{ExpressionElement, RuntimeAttributesSectionElement}
import wdl.model.draft3.graph.ExpressionValueConsumer
import wdl.model.draft3.graph.expression.{FileEvaluator, TypeEvaluator, ValueEvaluator}
import wdl.transforms.base.wdlom2wom.TaskDefinitionElementToWomTaskDefinition.{
  createTaskGraph,
  eliminateInputDependencies,
  expandLines,
  processMetaSections,
  validateParameterMetaEntries,
  TaskDefinitionElementToWomInputs
}
import wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.WorkflowDefinitionConvertInputs
import wdl.transforms.base.wdlom2wom.{
  CommandPartElementToWomCommandPart,
  FileElementToWomBundle,
  FileElementToWomBundleInputs,
  TaskDefinitionElementToWomTaskDefinition,
  WorkflowDefinitionElementToWomWorkflowDefinition
}
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wdl.transforms.cascades.linking.expression.consumed._
import wdl.transforms.cascades.linking.expression.files._
import wdl.transforms.cascades.linking.expression.types._
import wdl.transforms.cascades.linking.expression.values._
import wom.values.WomInteger
import wom.{CommandPart, RuntimeAttributes, RuntimeAttributesKeys}

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
    val inputElements = a.taskDefinitionElement.inputsSection.map(_.inputDeclarations).getOrElse(Seq.empty)

    val declarations = a.taskDefinitionElement.declarations
    val outputElements = a.taskDefinitionElement.outputsSection.map(_.outputs).getOrElse(Seq.empty)

    val conversion = (
      createTaskGraph(inputElements,
                      declarations,
                      outputElements,
                      a.taskDefinitionElement.parameterMetaSection,
                      a.typeAliases
      )(expressionValueConsumer, fileEvaluator, typeEvaluator, valueEvaluator),
      validateParameterMetaEntries(a.taskDefinitionElement.parameterMetaSection,
                                   a.taskDefinitionElement.inputsSection,
                                   a.taskDefinitionElement.outputsSection
      )
    ) flatMapN { (taskGraph, _) =>
      val validRuntimeAttributes: ErrorOr[RuntimeAttributes] = a.taskDefinitionElement.runtimeSection match {
        case Some(attributeSection) =>
          val returnCodesAttribute =
            attributeSection.runtimeAttributes.toList.find(pair =>
              pair.key.equals(RuntimeAttributesKeys.ReturnCodesKey)
            )
          val continueOnReturnCodeAttribute =
            attributeSection.runtimeAttributes.toList.find(pair =>
              pair.key.equals(RuntimeAttributesKeys.ContinueOnReturnCodeKey)
            )

          val returnCodesGet = returnCodesAttribute.orNull
          val continueOnReturnCodeGet = continueOnReturnCodeAttribute.orNull
          var editedAttributes = attributeSection.runtimeAttributes
          var returnCodesNotUnique = false

          if (returnCodesGet != null && continueOnReturnCodeGet != null) {
            returnCodesNotUnique = returnCodesGet.value
              .equals(continueOnReturnCodeGet.value) || returnCodesGet.value.equals(
              ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(0))))
            )
          }

          if (returnCodesGet != null && !returnCodesNotUnique) {
            editedAttributes = attributeSection.runtimeAttributes.filterNot(attribute =>
              attribute.key.equals(RuntimeAttributesKeys.ContinueOnReturnCodeKey)
            )
            editedAttributes = editedAttributes ++ Vector(
              KvPair(RuntimeAttributesKeys.ContinueOnReturnCodeKey, returnCodesAttribute.get.value)
            )
          }

          editedAttributes =
            editedAttributes.filterNot(attribute => attribute.key.equals(RuntimeAttributesKeys.ReturnCodesKey))

          TaskDefinitionElementToWomTaskDefinition.createRuntimeAttributes(
            RuntimeAttributesSectionElement(editedAttributes),
            taskGraph.linkedGraph
          )(
            expressionValueConsumer,
            fileEvaluator,
            typeEvaluator,
            valueEvaluator
          )
        case None => RuntimeAttributes(Map.empty).validNel
      }

      val validCommand: ErrorOr[Seq[CommandPart]] =
        expandLines(a.taskDefinitionElement.commandSection.parts).toList
          .traverse { parts =>
            CommandPartElementToWomCommandPart.convert(parts,
                                                       taskGraph.linkedGraph.typeAliases,
                                                       taskGraph.linkedGraph.generatedHandles
            )(expressionValueConsumer, fileEvaluator, typeEvaluator, valueEvaluator)
          }
          .map(_.toSeq)

      val (meta, parameterMeta) =
        processMetaSections(a.taskDefinitionElement.metaSection, a.taskDefinitionElement.parameterMetaSection)

      (validRuntimeAttributes, validCommand) mapN { (runtime, command) =>
        CallableTaskDefinition(
          a.taskDefinitionElement.name,
          Function.const(command.validNel),
          runtime,
          meta,
          parameterMeta,
          taskGraph.outputs,
          taskGraph.inputs,
          Set.empty,
          Map.empty,
          sourceLocation = a.taskDefinitionElement.sourceLocation
        )
      }
    }

    conversion.contextualizeErrors(s"process task definition '${b.taskDefinitionElement.name}'")
  }
}
