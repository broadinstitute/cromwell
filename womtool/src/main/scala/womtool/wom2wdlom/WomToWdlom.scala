package womtool.wom2wdlom

import cats.instances.either._
import cats.instances.list._
import cats.syntax.traverse._
import common.Checked
import common.collections.EnhancedCollections.EnhancedTraversableLike
import common.transforms.CheckedAtoB
import common.validation.Checked._
import shapeless.{Inl, Inr}
import wdl.draft2.model.command.{ParameterCommandPart, StringCommandPart}
import wdl.model.draft3.elements._
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement.ExpressionLiteralElement
import wom.callable.Callable._
import wom.callable.Callable.OutputDefinition
import wom.callable.{CallableTaskDefinition, MetaValueElement, WorkflowDefinition}
import wom.executable.WomBundle
import wom.expression.WomExpression
import wom.graph._
import wom.graph.CallNode.InputDefinitionPointer
import wom.graph.expression._
import wom.graph.GraphNodePort.GraphNodeOutputPort
import wom.RuntimeAttributes
import wom.types._

object WomToWdlom {

  private def invalidFromString(text: String) =
    s"Value \'$text\' has no representation in the destination format (WDL)".invalidNelCheck

  def graphOutputNodeToWorkflowGraphElement: CheckedAtoB[GraphOutputNode, WorkflowGraphElement] =
    CheckedAtoB.fromCheck { a: GraphOutputNode => a match {
      case a: ExpressionBasedGraphOutputNode =>
        for {
          typeElement <- womTypeToTypeElement(a.womType)
          expression <- womExpressionToExpressionElement(a.womExpression)
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

      for {
        workflows <- workflows.map(workflowDefinitionToWorkflowDefinitionElement(_)).toList.sequence[Checked, WorkflowDefinitionElement]
        tasks <- tasks.map(callableTaskDefinitionToTaskDefinitionElement(_)).toList.sequence[Checked, TaskDefinitionElement]
      } yield {
        FileElement(
          Seq.empty, // Imports do not exist in WOM and cannot sensibly be added at this point
          Seq.empty, // Structs do not exist in draft-2
          workflows,
          tasks
        )
      }
    }

  def mapToMetaSectionElement: CheckedAtoB[Map[String, MetaValueElement], Option[MetaSectionElement]] =
    CheckedAtoB.fromCheck { a: Map[String, MetaValueElement] =>
      if (a.nonEmpty)
        Some(MetaSectionElement(a)).validNelCheck
      else
        None.validNelCheck
    }

  def mapToParameterMetaSectionElement: CheckedAtoB[Map[String, MetaValueElement], Option[ParameterMetaSectionElement]] =
    CheckedAtoB.fromCheck { a: Map[String, MetaValueElement] =>
      if (a.nonEmpty)
        Some(ParameterMetaSectionElement(a)).validNelCheck
      else
        None.validNelCheck
    }

  def runtimeAttributesToRuntimeAttributesSectionElement: CheckedAtoB[RuntimeAttributes, Option[RuntimeAttributesSectionElement]] =
    CheckedAtoB.fromCheck { a: RuntimeAttributes =>
      def tupleToKvPair(tuple: (String, WomExpression)): Checked[ExpressionElement.KvPair] = {
        for {
          expressionElement <- womExpressionToExpressionElement(tuple._2)
        } yield {
          ExpressionElement.KvPair(tuple._1, expressionElement)
        }
      }

      for {
        kvPairs <- (a.attributes map tupleToKvPair).toList.sequence[Checked, ExpressionElement.KvPair]
      } yield {
        if (kvPairs.nonEmpty)
          Some(RuntimeAttributesSectionElement(kvPairs.toVector))
        else
          None
      }
    }

  def outputDefinitionToOutputDeclarationElement: CheckedAtoB[OutputDefinition, OutputDeclarationElement] =
    CheckedAtoB.fromCheck { a: OutputDefinition =>
      for {
        typeElement <- womTypeToTypeElement(a.womType)
        expression <- womExpressionToExpressionElement(a.expression)
      } yield {
        OutputDeclarationElement(typeElement, a.name, expression)
      }
    }

  def inputDefinitionToInputDeclarationElement: CheckedAtoB[InputDefinition, InputDeclarationElement] =
    CheckedAtoB.fromCheck { a: InputDefinition =>
      womTypeToTypeElement(a.womType) flatMap { womType =>
        a match {
          case a: RequiredInputDefinition =>
            InputDeclarationElement(womType, a.localName.value, None).validNelCheck
          case a: OverridableInputDefinitionWithDefault =>
            womExpressionToExpressionElement(a.default) map { expression =>
              InputDeclarationElement(womType, a.localName.value, Some(expression))
            }
          case a: FixedInputDefinitionWithDefault =>
            womExpressionToExpressionElement(a.default) map { expression =>
              InputDeclarationElement(womType, a.localName.value, Some(expression))
            }
          case a: OptionalInputDefinition =>
            InputDeclarationElement(womType, a.localName.value, None).validNelCheck
        }
      }
    }

  def callableTaskDefinitionToTaskDefinitionElement: CheckedAtoB[CallableTaskDefinition, TaskDefinitionElement] =
    CheckedAtoB.fromCheck {
      a: CallableTaskDefinition =>
        val inputs: Checked[List[InputDeclarationElement]] =
          a.inputs.map(inputDefinitionToInputDeclarationElement(_)).sequence[Checked, InputDeclarationElement]
        val outputs: Checked[List[OutputDeclarationElement]] =
          a.outputs.map(outputDefinitionToOutputDeclarationElement(_)).sequence[Checked, OutputDeclarationElement]

        for {
          runtime <- runtimeAttributesToRuntimeAttributesSectionElement(a.runtimeAttributes)
          meta <- mapToMetaSectionElement(a.meta)
          parameterMeta <- mapToParameterMetaSectionElement(a.parameterMeta)
          inputs <- inputs
          outputs <- outputs
          commands <- a.commandTemplateBuilder(Map()).toEither
        } yield {
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

          TaskDefinitionElement(
            a.name,
            if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
            Seq.empty, // No such thing in draft-2
            if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
            CommandSectionElement(Seq(commandLine)),
            runtime,
            meta,
            parameterMeta,
            a.sourceLocation
          )
        }
      }

  def workflowDefinitionToWorkflowDefinitionElement: CheckedAtoB[WorkflowDefinition, WorkflowDefinitionElement] =
    CheckedAtoB.fromCheck { a: WorkflowDefinition =>
        // This is a bit odd, so let's explain. "Real" inputs/outputs that are specified by the WDL's author
        // cannot have periods in them - period. So any input/output that has a period in it
        // is an artifact of WOMification and should be dropped
        val inputs =
          a.inputs.filter(!_.localName.value.contains(".")).map(inputDefinitionToInputDeclarationElement(_)).sequence[Checked, InputDeclarationElement]
        val outputs =
          a.outputs.filter(!_.localName.value.contains(".")).map(outputDefinitionToOutputDeclarationElement(_)).sequence[Checked, OutputDeclarationElement]

        for {
          meta <- mapToMetaSectionElement(a.meta)
          parameterMeta <- mapToParameterMetaSectionElement(a.parameterMeta)
          inputs <- inputs
          outputs <- outputs
          nodes <- selectWdlomRepresentableNodes(a.graph.nodes).map(graphNodeToWorkflowGraphElement(_)).toList.sequence[Checked, WorkflowGraphElement]
        } yield {
          WorkflowDefinitionElement(
            a.name,
            if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
            nodes.toSet,
            if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
            meta,
            parameterMeta,
            a.sourceLocation)
        }
      }

  def expressionNodeLikeToWorkflowGraphElement: CheckedAtoB[ExpressionNodeLike, WorkflowGraphElement] =
    CheckedAtoB.fromCheck {
      case a: ExpressionNode =>
        for {
          typeElement <- womTypeToTypeElement(a.womType)
          expression <- expressionNodeToExpressionElement(a)
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
        callNodeToCallElement(a)
      case a: ConditionalNode =>
        for {
          condition <- expressionNodeToExpressionElement(a.conditionExpression)
          nodes <- selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement(_)).sequence[Checked, WorkflowGraphElement]
        } yield {
          IfElement(
            conditionExpression = condition,
            graphElements = nodes
          )
        }
      case a: ExpressionNodeLike =>
        expressionNodeLikeToWorkflowGraphElement(a)
      case a: GraphNodeWithSingleOutputPort =>
        graphNodeWithSingleOutputPortToWorkflowGraphElement(a)
      case a: GraphOutputNode =>
        graphOutputNodeToWorkflowGraphElement(a)
      case a: ScatterNode =>
        // Only CWL has multi-variable scatters
        if (a.scatterVariableNodes.size != 1)
          s"In WDL, the scatter node count must be exactly one, was ${a.scatterVariableNodes.size}".invalidNelCheck
        else
          for {
            expression <- expressionNodeToExpressionElement(a.scatterCollectionExpressionNodes.head)
            graph <- selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement(_)).sequence[Checked, WorkflowGraphElement]
          } yield {
            ScatterElement(
              scatterName = a.identifier.localName.value,
              scatterExpression = expression,
              scatterVariableName = a.inputPorts.toList.head.name,
              graphElements = graph,
              sourceLocation = None
            )
          }
      }

  def graphNodeWithSingleOutputPortToWorkflowGraphElement: CheckedAtoB[GraphNodeWithSingleOutputPort, WorkflowGraphElement] =
    CheckedAtoB.fromCheck {
      case a: GraphInputNode =>
        womTypeToTypeElement(a.womType) map { typeElement =>
          InputDeclarationElement(
            typeElement,
            a.identifier.localName.value,
            None)
        }
      case a: ExpressionNode =>
        for {
          typeElement <- womTypeToTypeElement(a.womType)
          expression <- womExpressionToExpressionElement(a.womExpression)
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
        womTypeToTypeElement(a.memberType) map { typeElement =>
          if (a.guaranteedNonEmpty)
            NonEmptyTypeElement(ArrayTypeElement(typeElement))
          else
            ArrayTypeElement(typeElement)
        }
      case _: WomCoproductType => invalidFromString("WDL does not have coproducts - is this WOM from CWL?")
      case _: WomFileType => PrimitiveTypeElement(WomSingleFileType).validNelCheck
      case a: WomMapType =>
        for {
          keyType <- womTypeToTypeElement(a.keyType)
          valueType <- womTypeToTypeElement(a.valueType)
        } yield {
          MapTypeElement(keyType, valueType)
        }
      case _: WomNothingType.type => invalidFromString("WDL does not have the Nothing type - is this WOM from CWL?")
      case _: WomObjectType.type => ObjectTypeElement.validNelCheck
      case a: WomOptionalType => womOptionalTypeToOptionalTypeElement(a)
      case a: WomPairType =>
        for {
          leftType <- womTypeToTypeElement(a.leftType)
          rightType <- womTypeToTypeElement(a.rightType)
        } yield {
          PairTypeElement(leftType, rightType)
        }
      case a: WomPrimitiveType => womPrimitiveTypeToPrimitiveTypeElement(a)
    }

  def womOptionalTypeToOptionalTypeElement: CheckedAtoB[WomOptionalType, OptionalTypeElement] =
    CheckedAtoB.fromCheck { a: WomOptionalType =>
      womTypeToTypeElement(a.memberType) map { typeElement =>
        OptionalTypeElement(typeElement)
      }
    }

  def womPrimitiveTypeToPrimitiveTypeElement: CheckedAtoB[WomPrimitiveType, PrimitiveTypeElement] =
    CheckedAtoB.fromCheck { a: WomPrimitiveType =>
      PrimitiveTypeElement(a).validNelCheck
    }

  def expressionNodeToExpressionElement: CheckedAtoB[ExpressionNode, ExpressionElement] =
    CheckedAtoB.fromCheck {
      a: ExpressionNode => womExpressionToExpressionElement(a.womExpression)
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
          womExpressionToExpressionElement(a.womExpression) map { expressionElement =>
            Some(expressionElement)
          }
        case a: TaskCallInputExpressionNode =>
          womExpressionToExpressionElement(a.womExpression) map { expressionElement =>
            Some(expressionElement)
          }
        case _: RequiredGraphInputNode => None.validNelCheck
      }
      // Input definitions that directly contain expressions are the result of accepting a default input defined by the callable
      case Inr(Inl(_: WomExpression)) => None.validNelCheck
      case Inr(a) =>
        invalidFromString(a.toString)
      case a =>
        invalidFromString(a.toString)
    }

  def callNodeToCallElement: CheckedAtoB[CallNode, CallElement] =
    CheckedAtoB.fromCheck { call: CallNode =>
      def tupleToKvPair(tuple: (InputDefinition, InputDefinitionPointer)): Option[ExpressionElement.KvPair] = {
        inputDefinitionPointerToExpressionElement(tuple._2) match {
          case Right(Some(value)) => Some(ExpressionElement.KvPair(tuple._1.name, value))
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

      val afters = call.nonInputBasedPrerequisites.map(_.localName).toVector

      CallElement(
        callableName,
        maybeAlias,
        afters,
        if (inputs.nonEmpty) Some(CallBodyElement(inputs)) else None,
        call.sourceLocation
      ).validNelCheck
    }
}
