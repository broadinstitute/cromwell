package wdl.draft3.transforms.wdlom2wom

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.linking.graph.LinkedGraphMaker
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.LinkedGraph
import wom.{CommandPart, RuntimeAttributes}
import wom.callable.{Callable, CallableTaskDefinition}
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.draft3.transforms.linking.expression._
import wom.callable.Callable._
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.draft3.transforms.linking.typemakers._
import wdl.draft3.transforms.linking.expression.consumed._
import wdl.draft3.transforms.wdlom2wom.expression.renaming.IdentifierLookupRenamer.ops._
import wdl.draft3.transforms.wdlom2wom.expression.renaming.expressionEvaluator
import wdl.model.draft3.elements.ExpressionElement.{ArrayLiteral, IdentifierLookup, KvPair, SelectFirst}
import wom.expression.WomExpression
import wom.types.{WomOptionalType, WomType}

object TaskDefinitionElementToWomTaskDefinition {

  final case class TaskDefinitionElementToWomInputs(taskDefinitionElement: TaskDefinitionElement, typeAliases: Map[String, WomType])

  def convert(b: TaskDefinitionElementToWomInputs): ErrorOr[CallableTaskDefinition] = {
    val a = eliminateInputDependencies(b)
    val inputElements = a.taskDefinitionElement.inputsSection.map(_.inputDeclarations).getOrElse(Seq.empty)


    val declarations = a.taskDefinitionElement.declarations
    val outputElements = a.taskDefinitionElement.outputsSection.map(_.outputs).getOrElse(Seq.empty)

    val conversion = createTaskGraph(inputElements, declarations, outputElements, a.typeAliases) flatMap { taskGraph =>
      val validRuntimeAttributes: ErrorOr[RuntimeAttributes] = a.taskDefinitionElement.runtimeSection match {
        case Some(attributeSection) => createRuntimeAttributes(attributeSection, taskGraph.linkedGraph)
        case None => RuntimeAttributes(Map.empty).validNel
      }

      val validCommand: ErrorOr[Seq[CommandPart]] = {
        expandLines(a.taskDefinitionElement.commandSection.parts).toList.traverse { parts =>
          CommandPartElementToWomCommandPart.convert(parts, taskGraph.linkedGraph.typeAliases, taskGraph.linkedGraph.generatedHandles)
        }.map(_.toSeq)
      }

      (validRuntimeAttributes, validCommand) mapN { (runtime, command) =>
        CallableTaskDefinition(a.taskDefinitionElement.name, Function.const(command.validNel), runtime, Map.empty, Map.empty, taskGraph.outputs, taskGraph.inputs, Set.empty, Map.empty)
      }
    }

    conversion.contextualizeErrors(s"process task definition '${b.taskDefinitionElement.name}'")
  }

  private def eliminateInputDependencies(a: TaskDefinitionElementToWomInputs): TaskDefinitionElementToWomInputs = {
    case class NewInputElementsSet(original: InputDeclarationElement, newInput: InputDeclarationElement, newDeclaration: IntermediateValueDeclarationElement)

    val inputElementsWithUpstreams: Seq[NewInputElementsSet] = a.taskDefinitionElement.inputsSection.map(_.inputDeclarations).getOrElse(Seq.empty) collect {
      case ide @ InputDeclarationElement(typeElement,name, Some(expression)) if expression.expressionConsumedValueHooks.nonEmpty =>
        val input = InputDeclarationElement(OptionalTypeElement(typeElement), name, None)

        val selecterExpression = SelectFirst(ArrayLiteral(Seq(IdentifierLookup(name), expression)))
        val intermediate = IntermediateValueDeclarationElement(typeElement, s"__$name", selecterExpression)

        NewInputElementsSet(ide, input, intermediate)
    }

    if (inputElementsWithUpstreams.nonEmpty) {
      val newInputsSection: Option[InputsSectionElement] = a.taskDefinitionElement.inputsSection.map { inputsSection =>
        InputsSectionElement(inputsSection.inputDeclarations.filterNot(inputElementsWithUpstreams.map(_.original).contains) ++ inputElementsWithUpstreams.map(_.newInput))
      }
      val newDeclarationsSection: Seq[IntermediateValueDeclarationElement] = a.taskDefinitionElement.declarations ++ inputElementsWithUpstreams.map(_.newDeclaration)

      val newTaskDefinitionElement = {
        val withNewInputElements = a.taskDefinitionElement.copy(inputsSection = newInputsSection, declarations = newDeclarationsSection)
        val identifierRenames = inputElementsWithUpstreams.map(ie => ie.original.name -> ie.newDeclaration.name).toMap
        renameIdentifierAccesses(withNewInputElements, identifierRenames)
      }

      a.copy(taskDefinitionElement = newTaskDefinitionElement)

    } else a
  }

  private def renameIdentifierAccesses(taskDefinitionElement: TaskDefinitionElement, renames: Map[String, String]) = {

    val updatedInputs = taskDefinitionElement.inputsSection.map { inputsSection => InputsSectionElement(
      inputsSection.inputDeclarations map { id =>
        id.copy(expression = id.expression.map(_.renameIdentifiers(renames)))
      }
    )}

    val updatedOutputs = taskDefinitionElement.outputsSection.map { outputsSection => OutputsSectionElement(
      outputsSection.outputs map { od =>
        od.copy(expression = od.expression.renameIdentifiers(renames))
      }
    )}

    val updatedCommand = CommandSectionElement(
      taskDefinitionElement.commandSection.parts.map { line => CommandSectionLine(line.parts.map {
        case s: StringCommandPartElement => s
        case p: PlaceholderCommandPartElement => p.copy(expressionElement = p.expressionElement.renameIdentifiers(renames))
      })})

    val updatedRuntime = taskDefinitionElement.runtimeSection.map { runtimeSection => RuntimeAttributesSectionElement(
      runtimeSection.runtimeAttributes.map {
        case KvPair(key, value) => KvPair(key, value.renameIdentifiers(renames))
      }
    )}

    taskDefinitionElement.copy(
      inputsSection = updatedInputs,
      outputsSection = updatedOutputs,
      commandSection = updatedCommand,
      runtimeSection = updatedRuntime
    )
  }

  private def expandLines(lines: Seq[CommandSectionLine]): Seq[CommandPartElement] = {
    def expandNonFinalLine(line: CommandSectionLine): Seq[CommandPartElement] = if (line.parts.isEmpty) {
      Seq(StringCommandPartElement(System.lineSeparator))
    } else {
      val finalElements = line.parts.lastOption match {
        case Some(StringCommandPartElement(str)) => Seq(StringCommandPartElement(str + System.lineSeparator))
        case Some(other) => Seq(other, StringCommandPartElement(System.lineSeparator))
        case None => Seq(StringCommandPartElement(System.lineSeparator))
      }

      line.parts.init ++ finalElements
    }

    if (lines.isEmpty) {
      Seq(StringCommandPartElement(""))
    } else {
      (lines.init flatMap expandNonFinalLine) ++ lines.last.parts
    }
  }

  private final case class TaskGraph(inputs: List[Callable.InputDefinition], outputs: List[Callable.OutputDefinition], linkedGraph: LinkedGraph)

  private def createTaskGraph(inputs: Seq[InputDeclarationElement],
                              declarations: Seq[IntermediateValueDeclarationElement],
                              outputs: Seq[OutputDeclarationElement],
                              typeAliases: Map[String, WomType]): ErrorOr[TaskGraph] = {
    val combined: Set[WorkflowGraphElement] = (inputs ++ declarations ++ outputs).toSet
    LinkedGraphMaker.make(combined, Set.empty, typeAliases, Map.empty) flatMap { linked =>
      val ordered = LinkedGraphMaker.getOrdering(linked)

      def foldFunction(currentGraphValidation: ErrorOr[TaskGraph], next: WorkflowGraphElement): ErrorOr[TaskGraph] = {
        currentGraphValidation flatMap { accumulator => addToTaskGraph(next, accumulator) }
      }

      def addToTaskGraph(element: WorkflowGraphElement, accumulator: TaskGraph): ErrorOr[TaskGraph] = element match {
        case IntermediateValueDeclarationElement(womTypeElement, name, expression) =>
          val typeValidation = womTypeElement.determineWomType(linked.typeAliases)
          val expressionValidation = expression.makeWomExpression(linked.typeAliases, linked.consumedValueLookup)

          (typeValidation, expressionValidation) mapN { (womType, womExpression) =>
            accumulator.copy(inputs = accumulator.inputs :+ FixedInputDefinition(name, womType, womExpression))
          }
        case InputDeclarationElement(womTypeElement, name, None) =>
          womTypeElement.determineWomType(linked.typeAliases) map { womType =>
            val newInput = womType match {
              case optional: WomOptionalType => OptionalInputDefinition(name, optional)
              case required => RequiredInputDefinition(name, required)
            }
            accumulator.copy(inputs = accumulator.inputs :+ newInput)
          }
        case InputDeclarationElement(womTypeElement, name, Some(expression)) if expression.expressionConsumedValueHooks.isEmpty =>
          val typeValidation = womTypeElement.determineWomType(linked.typeAliases)
          val expressionValidation = expression.makeWomExpression(linked.typeAliases, linked.consumedValueLookup)

          (typeValidation, expressionValidation) mapN { (womType, womExpression) =>
            accumulator.copy(inputs = accumulator.inputs :+ InputDefinitionWithDefault(name, womType, womExpression))
          }

        // In this case, the expression has upstream dependencies. Since WOM won't allow that, make this an optional input and fixed declaration pair:
        case InputDeclarationElement(womTypeElement, name, Some(expression)) =>
          womTypeElement.determineWomType(linked.typeAliases) flatMap { womType =>
            val newInputName = s"__$name"
            val newInputType = WomOptionalType(womType).flatOptionalType
            val newInputDefinition = OptionalInputDefinition(newInputName, newInputType)

            val intermediateExpression: ExpressionElement = SelectFirst(ArrayLiteral(Seq(IdentifierLookup(newInputName), expression)))

            intermediateExpression.makeWomExpression(linked.typeAliases, linked.consumedValueLookup) map { womExpression =>
              val intermediateDefinition = FixedInputDefinition(name, womType, womExpression)

              accumulator.copy(inputs = accumulator.inputs :+ newInputDefinition :+ intermediateDefinition)
            }
          }


        case OutputDeclarationElement(womTypeElement, name, expression) =>
          val typeValidation = womTypeElement.determineWomType(linked.typeAliases)
          val expressionValidation = expression.makeWomExpression(linked.typeAliases, linked.consumedValueLookup)

          (typeValidation, expressionValidation) mapN { (womType, womExpression) =>
            accumulator.copy(outputs = accumulator.outputs :+ OutputDefinition(name, womType, womExpression))
          }
      }

      val initialState: ErrorOr[TaskGraph] = TaskGraph(List.empty, List.empty, linked).validNel
      ordered flatMap {
        _.foldLeft(initialState)(foldFunction)
      }
    }
  }

  private def createRuntimeAttributes(attributes: RuntimeAttributesSectionElement, linkedGraph: LinkedGraph): ErrorOr[RuntimeAttributes] = {

    def processSingleRuntimeAttribute(kvPair: KvPair): ErrorOr[(String, WomExpression)] = for {
      consumedValueLookup <- LinkedGraphMaker.makeConsumedValueLookup(kvPair.value.expressionConsumedValueHooks, linkedGraph.generatedHandles)
      womExpression <- kvPair.value.makeWomExpression(linkedGraph.typeAliases, consumedValueLookup)
    } yield kvPair.key -> womExpression


    attributes.runtimeAttributes.toList.traverse(processSingleRuntimeAttribute).map(atts => RuntimeAttributes(atts.toMap))
  }

}
