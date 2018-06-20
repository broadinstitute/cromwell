package womtool.wom2wdlom

import cats.data.Validated.Valid
import common.validation.Checked._
import wdl.model.draft3.elements._
import wom.executable.WomBundle
import common.collections.EnhancedCollections.EnhancedTraversableLike
import common.transforms.CheckedAtoB
import wdl.draft2.model.command.{ParameterCommandPart, StringCommandPart}
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import shapeless.{Inl, Inr}
import wom.callable.Callable._
import wdl.model.draft3.elements.ExpressionElement.ExpressionLiteralElement
import wom.callable.MetaValueElement.MetaValueElementString
import wom.RuntimeAttributes
import wom.callable.Callable.OutputDefinition
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.expression.WomExpression
import wom.graph.CallNode.InputDefinitionPointer
import wom.graph.GraphNodePort.GraphNodeOutputPort
import wom.graph._
import wom.graph.expression._
import wom.types._

case object UnrepresentableException extends Exception("Value has no representation in the destination format (WDL)")

object WomToWdlomImpl {

  private def invalidFromString(text: String) =
    s"Value \'$text\' has no representation in the destination format (WDL)".invalidNelCheck

  def graphOutputNodeToWorkflowGraphElement: CheckedAtoB[GraphOutputNode, WorkflowGraphElement] =
    CheckedAtoB.fromCheck { a: GraphOutputNode => a match {
      case a: ExpressionBasedGraphOutputNode =>
        for {
          typeElement <- womTypeToTypeElement.run(a.womType)
          expression <- womExpressionToExpressionElement.run(a.womExpression)
        } yield {
          OutputDeclarationElement(
            typeElement,
            a.identifier.localName.value,
            expression)
        }
      case _ =>
        invalidFromString(a.toString)
      }
    }

  def womBundleToFileElement: CheckedAtoB[WomBundle, FileElement] =
    CheckedAtoB.fromCheck { a: WomBundle =>
      val tasks: Iterable[CallableTaskDefinition] = a.allCallables.values.filterByType[CallableTaskDefinition]
      val workflows: Iterable[WorkflowDefinition] = a.allCallables.values.filterByType[WorkflowDefinition]

      FileElement(
        Seq.empty, // TODO: imports
        Seq.empty, // Structs do not exist in draft-2
        workflows.map(workflowDefinitionToWorkflowDefinitionElement.run(_).getOrElse(???)).toSeq,
        tasks.map(callableTaskDefinitionToTaskDefinitionElement.run(_).getOrElse(???)).toSeq
      ).validNelCheck
    }

  def mapToMetaSectionElement: CheckedAtoB[Map[String, String], Option[MetaSectionElement]] =
    CheckedAtoB.fromCheck { a: Map[String, String] =>
      if (a.nonEmpty)
        Some(MetaSectionElement(a map { case (key, value) =>
          key -> MetaValueElementString(value) // draft-2: strings only
        })).validNelCheck
      else
        None.validNelCheck
    }

  def mapToParameterMetaSectionElement: CheckedAtoB[Map[String, String], Option[ParameterMetaSectionElement]] =
    CheckedAtoB.fromCheck { a: Map[String, String] =>
      if (a.nonEmpty)
        Some(ParameterMetaSectionElement(a map { case (key, value) =>
          key -> MetaValueElementString(value) // draft-2: strings only
        })).validNelCheck
      else
        None.validNelCheck
    }

  def runtimeAttributesToRuntimeAttributesSectionElement: CheckedAtoB[RuntimeAttributes, Option[RuntimeAttributesSectionElement]] =
    CheckedAtoB.fromCheck { a: RuntimeAttributes =>
      def tupleToKvPair(tuple: (String, WomExpression)): ExpressionElement.KvPair =
        ExpressionElement.KvPair(tuple._1, womExpressionToExpressionElement.run(tuple._2).getOrElse(???))

      val kvPairs = (a.attributes map tupleToKvPair).toVector

      if (kvPairs.nonEmpty)
        Some(RuntimeAttributesSectionElement(kvPairs)).validNelCheck
      else
        None.validNelCheck
    }

  def outputDefinitionToOutputDeclarationElement: CheckedAtoB[OutputDefinition, OutputDeclarationElement] =
    CheckedAtoB.fromCheck { a: OutputDefinition =>
      OutputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.name, womExpressionToExpressionElement.run(a.expression).getOrElse(???)).validNelCheck
    }

  def inputDefinitionToInputDeclarationElement: CheckedAtoB[InputDefinition, InputDeclarationElement] =
    CheckedAtoB.fromCheck {
      case a: RequiredInputDefinition =>
        InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, None).validNelCheck
      case a: InputDefinitionWithDefault =>
        InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, Some(womExpressionToExpressionElement.run(a.default).getOrElse(???))).validNelCheck
      case a: FixedInputDefinition =>
        InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, Some(womExpressionToExpressionElement.run(a.default).getOrElse(???))).validNelCheck
      case a: OptionalInputDefinition =>
        InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, None).validNelCheck
    }

  def callableTaskDefinitionToTaskDefinitionElement: CheckedAtoB[CallableTaskDefinition, TaskDefinitionElement] =
    CheckedAtoB.fromCheck {
      a: CallableTaskDefinition =>
        val inputs = a.inputs.map(inputDefinitionToInputDeclarationElement.run(_).getOrElse(???))
        val outputs = a.outputs.map(outputDefinitionToOutputDeclarationElement.run(_).getOrElse(???))

        val commands = a.commandTemplateBuilder(Map()) match {
          case Valid(cmd) => cmd
          case _ => throw UnrepresentableException
        }

        val commandLine = CommandSectionLine(commands map {
          case s: StringCommandPart =>
            StringCommandPartElement(s.literal)
          case p: ParameterCommandPart =>
            val attrs = PlaceholderAttributeSet(
              defaultAttribute = p.attributes.get("default"),
              trueAttribute = p.attributes.get("true"),
              falseAttribute = p.attributes.get("false"),
              sepAttribute = p.attributes.get("sep")
            )

            PlaceholderCommandPartElement(ExpressionLiteralElement(p.expression.toWomString), attrs)
        })

        for {
          runtime <- runtimeAttributesToRuntimeAttributesSectionElement.run(a.runtimeAttributes)
          meta <- mapToMetaSectionElement.run(a.meta)
          parameterMeta <- mapToParameterMetaSectionElement.run(a.parameterMeta)
        } yield {
          TaskDefinitionElement(
            a.name,
            if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
            Seq.empty, // No such thing in draft-2
            if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
            CommandSectionElement(Seq(commandLine)),
            runtime,
            meta,
            parameterMeta
          )
        }
      }

  def workflowDefinitionToWorkflowDefinitionElement: CheckedAtoB[WorkflowDefinition, WorkflowDefinitionElement] =
    CheckedAtoB.fromCheck { a: WorkflowDefinition =>
        // This is a bit odd, so let's explain. "Real" inputs/outputs that are specified by the WDL's author
        // cannot have periods in them - period. So any input/output that has a period in it
        // is an artifact of WOMification and should be dropped
        val inputs = a.inputs.filter(!_.localName.value.contains(".")).map(inputDefinitionToInputDeclarationElement.run(_).getOrElse(???))
        val outputs = a.outputs.filter(!_.localName.value.contains(".")).map(outputDefinitionToOutputDeclarationElement.run(_).getOrElse(???))

        for {
          meta <- mapToMetaSectionElement.run(a.meta)
          parameterMeta <- mapToParameterMetaSectionElement.run(a.parameterMeta)
        } yield {
          WorkflowDefinitionElement(
            a.name,
            if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
            selectWdlomRepresentableNodes(a.graph.nodes).map(graphNodeToWorkflowGraphElement.run(_).getOrElse(???)),
            if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
            meta,
            parameterMeta)
        }
      }

  def expressionNodeLikeToWorkflowGraphElement: CheckedAtoB[ExpressionNodeLike, WorkflowGraphElement] =
    CheckedAtoB.fromCheck {
      case a: ExpressionNode =>
        for {
          typeElement <- womTypeToTypeElement.run(a.womType)
          expression <- expressionNodeToExpressionElement.run(a)
        } yield {
          IntermediateValueDeclarationElement(
            typeElement = typeElement,
            name = a.identifier.localName.value,
            expression = expression)
        }
      case a: ExpressionCallNode => invalidFromString(a.toString)
    }

  // Select only nodes that have a Wdlom representation (i.e. were not synthesized during WOMification)
  // WOM has some explicit representations that are implicit in WDL; they are necessary for execution,
  // but do not make sense (or are illegal) in WDL source.
  private def selectWdlomRepresentableNodes(allNodes: Set[GraphNode]): Set[GraphNode] = {
    val expressions: Set[GraphNode] = allNodes.filterByType[ExposedExpressionNode]
    val scatters: Set[GraphNode] = allNodes.filterByType[ScatterNode]
    val calls: Set[GraphNode] = allNodes.filterByType[CallNode]
    val conditionals: Set[GraphNode] = allNodes.filterByType[ConditionalNode]

    expressions ++ scatters ++ calls ++ conditionals
  }

  def graphNodeToWorkflowGraphElement: CheckedAtoB[GraphNode, WorkflowGraphElement] =
    CheckedAtoB.fromCheck {
        case a: CallNode =>
          callNodeToCallElement.run(a)
        case a: ConditionalNode =>
          IfElement(
            conditionExpression = expressionNodeToExpressionElement.run(a.conditionExpression).getOrElse(???),
            graphElements = selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement.run(_).getOrElse(???))
          ).validNelCheck
        case a: ExpressionNodeLike =>
          expressionNodeLikeToWorkflowGraphElement.run(a)
        case a: GraphNodeWithSingleOutputPort =>
          graphNodeWithSingleOutputPortToWorkflowGraphElement.run(a)
        case a: GraphOutputNode =>
          graphOutputNodeToWorkflowGraphElement.run(a)
        case a: ScatterNode =>
          // Only CWL has multi-variable scatters
          if (a.scatterVariableNodes.size != 1) throw UnrepresentableException // TODO: upgrade from exceptions to typed errors

          ScatterElement(
            scatterName = a.identifier.localName.value,
            scatterExpression = expressionNodeToExpressionElement.run(a.scatterCollectionExpressionNodes.head).getOrElse(???),
            scatterVariableName = a.inputPorts.toList.head.name,
            graphElements = selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement.run(_).getOrElse(???))
          ).validNelCheck
      }

  def graphNodeWithSingleOutputPortToWorkflowGraphElement: CheckedAtoB[GraphNodeWithSingleOutputPort, WorkflowGraphElement] =
    CheckedAtoB.fromCheck {
      case a: GraphInputNode =>
        womTypeToTypeElement.run(a.womType) map { typeElement =>
          InputDeclarationElement(
            typeElement,
            a.identifier.localName.value,
            None)
        }
      case a: ExpressionNode =>
        for {
          typeElement <- womTypeToTypeElement.run(a.womType)
          expression <- womExpressionToExpressionElement.run(a.womExpression)
        } yield {
          IntermediateValueDeclarationElement(
            typeElement,
            a.identifier.localName.value,
            expression)
        }
    }

  def womTypeToTypeElement: CheckedAtoB[WomType, TypeElement] =
    CheckedAtoB.fromCheck {
      case a: WomArrayType =>
        womTypeToTypeElement.run(a.memberType) map { typeElement =>
          if (a.guaranteedNonEmpty)
            NonEmptyTypeElement(ArrayTypeElement(typeElement))
          else
            ArrayTypeElement(typeElement)
        }
      case _: WomCoproductType => throw UnrepresentableException
      case _: WomFileType => PrimitiveTypeElement(WomSingleFileType).validNelCheck
      case a: WomMapType =>
        for {
          keyType <- womTypeToTypeElement.run(a.keyType)
          valueType <- womTypeToTypeElement.run(a.valueType)
        } yield {
          MapTypeElement(keyType, valueType)
        }
      case _: WomNothingType.type => throw UnrepresentableException
      case _: WomObjectType.type => ObjectTypeElement.validNelCheck
      case a: WomOptionalType => womOptionalTypeToOptionalTypeElement.run(a)
      case a: WomPairType =>
        for {
          leftType <- womTypeToTypeElement.run(a.leftType)
          rightType <- womTypeToTypeElement.run(a.rightType)
        } yield {
          PairTypeElement(leftType, rightType)
        }
      case a: WomPrimitiveType => womPrimitiveTypeToPrimitiveTypeElement.run(a)
    }

  def womOptionalTypeToOptionalTypeElement: CheckedAtoB[WomOptionalType, OptionalTypeElement] =
    CheckedAtoB.fromCheck { a: WomOptionalType =>
      womTypeToTypeElement.run(a.memberType) map { typeElement =>
        OptionalTypeElement(typeElement)
      }
    }

  def womPrimitiveTypeToPrimitiveTypeElement: CheckedAtoB[WomPrimitiveType, PrimitiveTypeElement] =
    CheckedAtoB.fromCheck { a: WomPrimitiveType =>
      PrimitiveTypeElement(a).validNelCheck
    }

  def expressionNodeToExpressionElement: CheckedAtoB[ExpressionNode, ExpressionElement] =
    CheckedAtoB.fromCheck {
      a: ExpressionNode => womExpressionToExpressionElement.run(a.womExpression)
    }

  def womExpressionToExpressionElement: CheckedAtoB[WomExpression, ExpressionElement] =
    CheckedAtoB.fromCheck {
      case a: WomExpression => ExpressionLiteralElement(a.sourceString).validNelCheck
      case a => invalidFromString(a.toString)
    }

  def inputDefinitionPointerToExpressionElement: CheckedAtoB[InputDefinitionPointer, Option[ExpressionElement]] =
    CheckedAtoB.fromCheck {
      // If the input definition is a node containing an expression, it's been declared explicitly
      case Inl(a: GraphNodeOutputPort) => a.graphNode match {
        case _: OptionalGraphInputNode =>
          None.validNelCheck
        case _: OptionalGraphInputNodeWithDefault =>
          None.validNelCheck
        case a: PlainAnonymousExpressionNode =>
          womExpressionToExpressionElement.run(a.womExpression) map { expressionElement =>
            Some(expressionElement)
          }
        case a: TaskCallInputExpressionNode =>
          womExpressionToExpressionElement.run(a.womExpression) map { expressionElement =>
            Some(expressionElement)
          }
        case _: RequiredGraphInputNode => None.validNelCheck
      }
      // Input definitions that directly contain expressions are the result of accepting a default input defined by the callable
      case Inr(Inl(_: WomExpression)) => None.validNelCheck
      case Inr(_) =>
        throw UnrepresentableException
      case a =>
        invalidFromString(a.toString)
    }

  def callNodeToCallElement: CheckedAtoB[CallNode, CallElement] =
    CheckedAtoB.fromCheck { call: CallNode =>
      def tupleToKvPair(tuple: (InputDefinition, InputDefinitionPointer)): Option[ExpressionElement.KvPair] = {
        inputDefinitionPointerToExpressionElement.run(tuple._2).getOrElse(???) match {
          case Some(value) => Some(ExpressionElement.KvPair(tuple._1.name, value))
          case _ => None
        }
      }

      val inputs = (call.inputDefinitionMappings flatMap tupleToKvPair).toVector

      val callableName = call.callable.name
      val callAlias = call.identifier.localName.value // If no alias, this is just the name; we evaluate for that below

      val maybeAlias = if (callableName != callAlias)
        Some(callAlias)
      else
        None

      CallElement(
        call.callable.name,
        maybeAlias,
        if (inputs.nonEmpty) Some(CallBodyElement(inputs)) else None
      ).validNelCheck
    }
}
