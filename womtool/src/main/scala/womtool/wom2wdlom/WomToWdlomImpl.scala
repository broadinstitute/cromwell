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
        OutputDeclarationElement(
          a.womType.toWdlom,
          a.identifier.localName.value,
          a.womExpression.toWdlom).validNelCheck
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

  implicit val mapToParameterMetaSectionElement: WomToWdlom[Map[String, String], Option[ParameterMetaSectionElement]] =
    new WomToWdlom[Map[String, String], Option[ParameterMetaSectionElement]] {
      override def toWdlom(a: Map[String, String]): Option[ParameterMetaSectionElement] = {
        if (a.nonEmpty)
          Some(ParameterMetaSectionElement(a map { case (key, value) =>
            key -> MetaValueElementString(value) // draft-2: strings only
          }))
        else
          None
      }
    }

  def runtimeAttributesToRuntimeAttributesSectionElement: CheckedAtoB[RuntimeAttributes, Option[RuntimeAttributesSectionElement]] =
    CheckedAtoB.fromCheck { a: RuntimeAttributes =>
      def tupleToKvPair(tuple: (String, WomExpression)): ExpressionElement.KvPair =
        ExpressionElement.KvPair(tuple._1, tuple._2.toWdlom)

      val kvPairs = (a.attributes map tupleToKvPair).toVector

      if (kvPairs.nonEmpty)
        Some(RuntimeAttributesSectionElement(kvPairs)).validNelCheck
      else
        None.validNelCheck
    }

  implicit val outputDefinitionToOutputDeclarationElement: WomToWdlom[OutputDefinition, OutputDeclarationElement] =
    new WomToWdlom[OutputDefinition, OutputDeclarationElement] {
      override def toWdlom(a: OutputDefinition): OutputDeclarationElement =
        OutputDeclarationElement(a.womType.toWdlom, a.name, a.expression.toWdlom)
    }

  implicit val inputDefinitionToInputDeclarationElement: WomToWdlom[InputDefinition, InputDeclarationElement] =
    new WomToWdlom[InputDefinition, InputDeclarationElement] {
      override def toWdlom(a: InputDefinition): InputDeclarationElement = a match {
        case a: RequiredInputDefinition => InputDeclarationElement(a.womType.toWdlom, a.localName.value, None)
        case a: InputDefinitionWithDefault => InputDeclarationElement(a.womType.toWdlom, a.localName.value, Some(a.default.toWdlom))
        case a: FixedInputDefinition => InputDeclarationElement(a.womType.toWdlom, a.localName.value, Some(a.default.toWdlom))
        case a: OptionalInputDefinition => InputDeclarationElement(a.womType.toWdlom, a.localName.value, None)
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
          mapToParameterMetaSectionElement.toWdlom(a.parameterMeta)
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
          mapToParameterMetaSectionElement.toWdlom(a.parameterMeta)
        )
      }
    }

  def expressionNodeLikeToWorkflowGraphElement: CheckedAtoB[ExpressionNodeLike, WorkflowGraphElement] =
    CheckedAtoB.fromCheck {
      case a: ExpressionNode =>
        IntermediateValueDeclarationElement(
          typeElement = a.womType.toWdlom,
          name = a.identifier.localName.value,
          expression = a.toWdlom).validNelCheck
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
            conditionExpression = a.conditionExpression.toWdlom,
            graphElements = selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement.toWdlom)
          )
        case a: ExpressionNodeLike =>
          expressionNodeLikeToWorkflowGraphElement.run(a).getOrElse(???)
        case a: GraphNodeWithSingleOutputPort =>
          a.toWdlom
        case a: GraphOutputNode =>
          graphOutputNodeToWorkflowGraphElement.run(a).getOrElse(???) // TODO: as de-exceptioning progresses, return bare Checked[]
        case a: ScatterNode =>
          // Only CWL has multi-variable scatters
          if (a.scatterVariableNodes.size != 1) throw UnrepresentableException // TODO: upgrade from exceptions to typed errors

          ScatterElement(
            scatterName = a.identifier.localName.value,
            scatterExpression = a.scatterCollectionExpressionNodes.head.toWdlom,
            scatterVariableName = a.inputPorts.toList.head.name,
            graphElements = selectWdlomRepresentableNodes(a.innerGraph.nodes).toList.map(graphNodeToWorkflowGraphElement.toWdlom)
          )
      }
    }
  }

  implicit val graphNodeWithSingleOutputPortToWorkflowGraphElement: WomToWdlom[GraphNodeWithSingleOutputPort, WorkflowGraphElement] =
    new WomToWdlom[GraphNodeWithSingleOutputPort, WorkflowGraphElement] {
      override def toWdlom(a: GraphNodeWithSingleOutputPort): WorkflowGraphElement = a match {
        case a: GraphInputNode =>
          InputDeclarationElement(
            a.womType.toWdlom,
            a.identifier.localName.value,
            None
          )
        case a: ExpressionNode =>
          IntermediateValueDeclarationElement(
            a.womType.toWdlom,
            a.identifier.localName.value,
            a.womExpression.toWdlom
          )
      }
    }

  implicit val womTypeToTypeElement: WomToWdlom[WomType, TypeElement] = new WomToWdlom[WomType, TypeElement] {
    override def toWdlom(a: WomType): TypeElement = a match {
      case a: WomArrayType =>
        if (a.guaranteedNonEmpty)
          NonEmptyTypeElement(ArrayTypeElement(womTypeToTypeElement.toWdlom(a.memberType)))
        else
          ArrayTypeElement(womTypeToTypeElement.toWdlom(a.memberType))
      case _: WomCoproductType => throw UnrepresentableException
      case _: WomFileType => PrimitiveTypeElement(WomSingleFileType)
      case a: WomMapType => MapTypeElement(womTypeToTypeElement.toWdlom(a.keyType), womTypeToTypeElement.toWdlom(a.valueType))
      case _: WomNothingType.type => throw UnrepresentableException
      case _: WomObjectType.type => ObjectTypeElement
      case a: WomOptionalType => a.toWdlom
      case a: WomPairType => PairTypeElement(womTypeToTypeElement.toWdlom(a.leftType), womTypeToTypeElement.toWdlom(a.rightType))
      case a: WomPrimitiveType => a.toWdlom
    }
  }

  implicit val womOptionalTypeToOptionalTypeElement: WomToWdlom[WomOptionalType, OptionalTypeElement] = new WomToWdlom[WomOptionalType, OptionalTypeElement] {
    override def toWdlom(a: WomOptionalType): OptionalTypeElement = OptionalTypeElement(womTypeToTypeElement.toWdlom(a.memberType))
  }

  implicit val womPrimitiveTypeToPrimitiveTypeElement: WomToWdlom[WomPrimitiveType, PrimitiveTypeElement] = new WomToWdlom[WomPrimitiveType, PrimitiveTypeElement] {
    override def toWdlom(a: WomPrimitiveType): PrimitiveTypeElement = PrimitiveTypeElement(a)
  }

  implicit val expressionNodeToExpressionElement: WomToWdlom[ExpressionNode, ExpressionElement] = new WomToWdlom[ExpressionNode, ExpressionElement] {
    override def toWdlom(a: ExpressionNode): ExpressionElement = a.womExpression.toWdlom
  }

  implicit val womExpressionToExpressionElement: WomToWdlom[WomExpression, ExpressionElement] = new WomToWdlom[WomExpression, ExpressionElement] {
    override def toWdlom(a: WomExpression): ExpressionElement = a match {
      case a: WomExpression => ExpressionLiteralElement(a.sourceString)
      case _ => throw UnrepresentableException
    }
  }

  def inputDefinitionPointerToExpressionElement: CheckedAtoB[InputDefinitionPointer, Option[ExpressionElement]] =
    CheckedAtoB.fromCheck {
      // If the input definition is a node containing an expression, it's been declared explicitly
      case Inl(a: GraphNodeOutputPort) => a.graphNode match {
        case _: OptionalGraphInputNode =>
          None.validNelCheck
        case _: OptionalGraphInputNodeWithDefault =>
          None.validNelCheck
        case a: PlainAnonymousExpressionNode => Some(a.womExpression.toWdlom).validNelCheck
        case a: TaskCallInputExpressionNode => Some(a.womExpression.toWdlom).validNelCheck
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
