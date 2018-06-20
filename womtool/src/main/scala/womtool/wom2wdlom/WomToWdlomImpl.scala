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

  import WomToWdlom.ops._

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
        workflows.map(_.toWdlom).toSeq,
        tasks.map(_.toWdlom).toSeq
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

  implicit val outputDefinitionToOutputDeclarationElement: WomToWdlom[OutputDefinition, OutputDeclarationElement] =
    new WomToWdlom[OutputDefinition, OutputDeclarationElement] {
      override def toWdlom(a: OutputDefinition): OutputDeclarationElement =
        OutputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.name, womExpressionToExpressionElement.run(a.expression).getOrElse(???))
    }

  implicit val inputDefinitionToInputDeclarationElement: WomToWdlom[InputDefinition, InputDeclarationElement] =
    new WomToWdlom[InputDefinition, InputDeclarationElement] {
      override def toWdlom(a: InputDefinition): InputDeclarationElement = a match {
        case a: RequiredInputDefinition =>
          InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, None)
        case a: InputDefinitionWithDefault =>
          InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, Some(womExpressionToExpressionElement.run(a.default).getOrElse(???)))
        case a: FixedInputDefinition =>
          InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, Some(womExpressionToExpressionElement.run(a.default).getOrElse(???)))
        case a: OptionalInputDefinition =>
          InputDeclarationElement(womTypeToTypeElement.run(a.womType).getOrElse(???), a.localName.value, None)
      }
    }

  implicit val callableTaskDefinitionToTaskDefinitionElement: WomToWdlom[CallableTaskDefinition, TaskDefinitionElement] =
    new WomToWdlom[CallableTaskDefinition, TaskDefinitionElement] {
      override def toWdlom(a: CallableTaskDefinition): TaskDefinitionElement = {
        val inputs = a.inputs.map(inputDefinitionToInputDeclarationElement.toWdlom)
        val outputs = a.outputs.map(outputDefinitionToOutputDeclarationElement.toWdlom)

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

        TaskDefinitionElement(
          a.name,
          if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
          Seq.empty, // No such thing in draft-2
          if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
          CommandSectionElement(Seq(commandLine)),
          runtimeAttributesToRuntimeAttributesSectionElement.run(a.runtimeAttributes).getOrElse(???),
          mapToMetaSectionElement.run(a.meta).getOrElse(???),
          mapToParameterMetaSectionElement.run(a.parameterMeta).getOrElse(???)
        )
      }
    }

  implicit val workflowDefinitionToWorkflowDefinitionElement: WomToWdlom[WorkflowDefinition, WorkflowDefinitionElement] =
    new WomToWdlom[WorkflowDefinition, WorkflowDefinitionElement] {
      override def toWdlom(a: WorkflowDefinition): WorkflowDefinitionElement = {
        // This is a bit odd, so let's explain. "Real" inputs/outputs that are specified by the WDL's author
        // cannot have periods in them - period. So any input/output that has a period in it
        // is an artifact of WOMification and should be dropped
        val inputs = a.inputs.filter(!_.localName.value.contains(".")).map(inputDefinitionToInputDeclarationElement.toWdlom)
        val outputs = a.outputs.filter(!_.localName.value.contains(".")).map(outputDefinitionToOutputDeclarationElement.toWdlom)

        WorkflowDefinitionElement(
          a.name,
          if (inputs.nonEmpty) Some(InputsSectionElement(inputs)) else None,
          selectWdlomRepresentableNodes(a.graph.nodes).map(_.toWdlom),
          if (outputs.nonEmpty) Some(OutputsSectionElement(outputs)) else None,
          mapToMetaSectionElement.run(a.meta).getOrElse(???),
          mapToParameterMetaSectionElement.run(a.parameterMeta).getOrElse(???)
        )
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


  implicit val graphNodeToWorkflowGraphElement: WomToWdlom[GraphNode, WorkflowGraphElement] = new WomToWdlom[GraphNode, WorkflowGraphElement] {
    override def toWdlom(a: GraphNode): WorkflowGraphElement = {
      a match {
        case a: CallNode =>
          callNodeToCallElement.run(a).getOrElse(???)
        case a: ConditionalNode =>
          IfElement(
            conditionExpression = expressionNodeToExpressionElement.run(a.conditionExpression).getOrElse(???),
            graphElements = selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement.toWdlom)
          )
        case a: ExpressionNodeLike =>
          expressionNodeLikeToWorkflowGraphElement.run(a).getOrElse(???)
        case a: GraphNodeWithSingleOutputPort =>
          graphNodeWithSingleOutputPortToWorkflowGraphElement.run(a).getOrElse(???)
        case a: GraphOutputNode =>
          graphOutputNodeToWorkflowGraphElement.run(a).getOrElse(???) // TODO: as de-exceptioning progresses, return bare Checked[]
        case a: ScatterNode =>
          // Only CWL has multi-variable scatters
          if (a.scatterVariableNodes.size != 1) throw UnrepresentableException // TODO: upgrade from exceptions to typed errors

          ScatterElement(
            scatterName = a.identifier.localName.value,
            scatterExpression = expressionNodeToExpressionElement.run(a.scatterCollectionExpressionNodes.head).getOrElse(???),
            scatterVariableName = a.inputPorts.toList.head.name,
            graphElements = selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement.toWdlom)
          )
      }
    }
  }

  def graphNodeWithSingleOutputPortToWorkflowGraphElement: CheckedAtoB[GraphNodeWithSingleOutputPort, WorkflowGraphElement] =
    CheckedAtoB.fromCheck {
      case a: GraphInputNode =>
        InputDeclarationElement(
          womTypeToTypeElement.run(a.womType).getOrElse(???),
          a.identifier.localName.value,
          None
        ).validNelCheck
      case a: ExpressionNode =>
        IntermediateValueDeclarationElement(
          womTypeToTypeElement.run(a.womType).getOrElse(???),
          a.identifier.localName.value,
          womExpressionToExpressionElement.run(a.womExpression).getOrElse(???)
        ).validNelCheck
    }

  def womTypeToTypeElement: CheckedAtoB[WomType, TypeElement] =
    CheckedAtoB.fromCheck {
      case a: WomArrayType =>
        if (a.guaranteedNonEmpty)
          NonEmptyTypeElement(ArrayTypeElement(womTypeToTypeElement.run(a.memberType).getOrElse(???))).validNelCheck
        else
          ArrayTypeElement(womTypeToTypeElement.run(a.memberType).getOrElse(???)).validNelCheck
      case _: WomCoproductType => throw UnrepresentableException
      case _: WomFileType => PrimitiveTypeElement(WomSingleFileType).validNelCheck
      case a: WomMapType => MapTypeElement(womTypeToTypeElement.run(a.keyType).getOrElse(???), womTypeToTypeElement.run(a.valueType).getOrElse(???)).validNelCheck
      case _: WomNothingType.type => throw UnrepresentableException
      case _: WomObjectType.type => ObjectTypeElement.validNelCheck
      case a: WomOptionalType => womOptionalTypeToOptionalTypeElement.run(a)
      case a: WomPairType => PairTypeElement(womTypeToTypeElement.run(a.leftType).getOrElse(???), womTypeToTypeElement.run(a.rightType).getOrElse(???)).validNelCheck
      case a: WomPrimitiveType => womPrimitiveTypeToPrimitiveTypeElement.run(a)
    }

  def womOptionalTypeToOptionalTypeElement: CheckedAtoB[WomOptionalType, OptionalTypeElement] =
    CheckedAtoB.fromCheck { a: WomOptionalType =>
      OptionalTypeElement(womTypeToTypeElement.run(a.memberType).getOrElse(???)).validNelCheck
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
        case a: PlainAnonymousExpressionNode => Some(womExpressionToExpressionElement.run(a.womExpression).getOrElse(???)).validNelCheck
        case a: TaskCallInputExpressionNode => Some(womExpressionToExpressionElement.run(a.womExpression).getOrElse(???)).validNelCheck
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
