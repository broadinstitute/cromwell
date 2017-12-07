package cwl

import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cats.syntax.validated._
import cwl.WorkflowStep.{ScatterMappings, WorkflowStepInputFold}
import shapeless.{Coproduct, Poly1}
import wom.callable.Callable.InputDefinition
import wom.graph.CallNode.{CallNodeAndNewNodes, CallNodeBuilder, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.types.WomArrayType

object ScatterLogic {

  object ScatterVariablesPoly extends Poly1 {
    implicit def fromString: Case.Aux[String, List[String]] = at[String] { s: String => List(s) }
    implicit def fromStringList: Case.Aux[Array[String], List[String]] = at[Array[String]] { l: Array[String] => l.toList }
  }

  // Validate that the scatter expression is an ArrayLike and return the member type
  def scatterExpressionItemType(expressionNode: ExpressionNode) = {
    expressionNode.womType match {
      case WomArrayType(itemType) => itemType.validNelCheck // Covers maps because this is a custom unapply (see WdlArrayType)
      case other => s"Cannot scatter over a non-traversable type ${other.toDisplayString}".invalidNelCheck
    }
  }

  // Build a list (potentially empty) of scatter helpers containing the expression node and scatter variable node
  // for each input being scattered over
  def buildScatterMappings(scatter: ScatterVariables, stepInputFold: WorkflowStepInputFold, stepId: String): Checked[ScatterMappings] = {
    import cats.implicits._

    def buildScatterHelperFor(scatterVariableName: String): Checked[(ExpressionNode, ScatterVariableNode)] = {
      // Assume the variable is a step input (is that always true ??). Find the corresponding expression node
      val parsedScatterVariable = FileStepAndId(scatterVariableName)

      stepInputFold.stepInputMapping.get(parsedScatterVariable.id) match {
        case Some(expressionNode) =>
          // find the item type
          scatterExpressionItemType(expressionNode) map { itemType =>
            // create a scatter variable node for other scattered nodes to point to
            val scatterVariableNode = ScatterVariableNode(WomIdentifier(scatterVariableName), expressionNode.singleExpressionOutputPort, itemType)
            (expressionNode, scatterVariableNode)
          }
        case None => s"Could not find a variable ${parsedScatterVariable.id} in the workflow step input to scatter over. Please make sure ${parsedScatterVariable.id} is an input of step $stepId".invalidNelCheck
      }
    }
    // Take the scatter field defining the (list of) input(s) to be scattered over
    scatter
      // Fold them so we have a List[String] no matter what (a single scattered input becomes a single element list)
      .map(_.fold(ScatterVariablesPoly))
      // If there's no scatter make it an empty list
      .getOrElse(List.empty)
      // Traverse the list to create ScatterHelpers that will be used later on to create the ScatterNode
      .traverse[Checked, (ExpressionNode, ScatterVariableNode)](buildScatterHelperFor)
      // Transform the List of tuples to a Map
      .map(_.toMap)
  }

  def buildScatteredInputFold(inputDefinition: InputDefinition,
                             scatterMappings: ScatterMappings,
                             expressionNode: ExpressionNode,
                             callNodeBuilder: CallNodeBuilder): ErrorOr[InputDefinitionFold] = {
    val (mapping, outputPort, newNode) = {
      // mapping = Where does this input definition get its value from (in this method it'll always be an output port, but wrapped in a Coproduct[InputDefinitionPointer]
      // outputPort = output port from which the input definition will get its values
      // newNode = potentially newly created OGIN to be added to the fold
      scatterMappings.get(expressionNode) match {
        // This input is being scattered over
        case Some(scatterVariableNode) =>
          // Point the input definition to the scatter variable
          val mapping = List(inputDefinition -> Coproduct[InputDefinitionPointer](scatterVariableNode.singleOutputPort: OutputPort))
          val outputPort = scatterVariableNode.singleOutputPort
          (mapping, outputPort, None)
        // This input is NOT being scattered over
        case None =>
          /*
            * If this input is not being scattered, we need to point at the step input expression node.
            * However the call node will be in the scatter inner graph, whereas the input expression node is outside of it.
            * So we create an OuterGraphInputNode to link them together.
           */
          val ogin = OuterGraphInputNode(WomIdentifier(inputDefinition.name), expressionNode.singleExpressionOutputPort, preserveScatterIndex = false)
          val mapping = List(inputDefinition -> Coproduct[InputDefinitionPointer](ogin.singleOutputPort: OutputPort))

          (mapping, ogin.singleOutputPort, Option(ogin))
      }
    }

    InputDefinitionFold(
      mappings = mapping,
      callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, outputPort)),
      newExpressionNodes = Set(expressionNode),
      usedOuterGraphInputNodes = newNode.toSet
    ).validNel
  }

  // Prepare the nodes to be returned if this call is being scattered
  def prepareNodesForScatteredCall(knownNodes: Set[GraphNode], callNodeAndNewNodes: CallNodeAndNewNodes, stepInputFold: WorkflowStepInputFold, scatterMappings: ScatterMappings): Checked[Set[GraphNode]] = {
    val callNode = callNodeAndNewNodes.node

    // We need to generate PBGONs for every output port of the call, so that they can be linked outside the scatter graph
    val portBasedGraphOutputNodes = callNode.outputPorts.map(op => {
      PortBasedGraphOutputNode(op.identifier, op.womType, op)
    })

    // Build the scatter inner graph using the callNode, the outerGraphInputNodes, the scatterVariableNodes, and the PBGONs
    Graph.validateAndConstruct(Set(callNodeAndNewNodes.node) ++ callNodeAndNewNodes.usedOuterGraphInputNodes ++ scatterMappings.values ++ portBasedGraphOutputNodes).toEither map { innerGraph =>
      // Build the ScatterNode - only with single scatter variable for now
      val scatterNodeWithNewNodes = ScatterNode.scatterOverGraph(innerGraph, scatterMappings.head._1, scatterMappings.head._2)
      knownNodes ++ scatterNodeWithNewNodes.nodes ++ stepInputFold.generatedNodes
    }
  }
}
