package cwl

import cats.data.NonEmptyList
import cats.syntax.validated._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cwl.ScatterMethod.ScatterMethod
import cwl.WorkflowStep.WorkflowStepInputFold
import shapeless.{Coproduct, Poly1}
import wom.callable.Callable.InputDefinition
import wom.graph.CallNode.{CallNodeAndNewNodes, CallNodeBuilder, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.ScatterNode.{ScatterCollectionFunctionBuilder, ScatterProcessingFunction, ScatterVariableAndValue}
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.types.WomArrayType
import wom.values.{WomArray, WomValue}

import scala.annotation.tailrec

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
  def buildScatterVariables(scatter: ScatterVariables, stepInputFold: WorkflowStepInputFold, stepId: String): Checked[List[ScatterVariableNode]] = {
    import cats.implicits._

    def buildScatterVariable(scatterVariableName: String): Checked[ScatterVariableNode] = {
      // Assume the variable is a step input (is that always true ??). Find the corresponding expression node
      val parsedScatterVariable = FileStepAndId(scatterVariableName)

      stepInputFold.stepInputMapping.get(parsedScatterVariable.id) match {
        case Some(expressionNode) =>
          // find the item type
          scatterExpressionItemType(expressionNode) map { itemType =>
            // create a scatter variable node for other scattered nodes to point to
            ScatterVariableNode(WomIdentifier(scatterVariableName), expressionNode, itemType)
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
      .traverse[Checked, ScatterVariableNode](buildScatterVariable)
  }

  def buildScatteredInputFold(inputDefinition: InputDefinition,
                              scatterMappings: List[ScatterVariableNode],
                              expressionNode: ExpressionNode,
                              callNodeBuilder: CallNodeBuilder): ErrorOr[InputDefinitionFold] = {
    val (mapping, outputPort, newNode) = {
      // mapping = Where does this input definition get its value from (in this method it'll always be an output port, but wrapped in a Coproduct[InputDefinitionPointer]
      // outputPort = output port from which the input definition will get its values
      // newNode = potentially newly created OGIN to be added to the fold
      scatterMappings.find(_.scatterExpressionNode == expressionNode) match {
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
  def prepareNodesForScatteredCall(knownNodes: Set[GraphNode],
                                   callNodeAndNewNodes: CallNodeAndNewNodes,
                                   stepInputFold: WorkflowStepInputFold,
                                   scatterVariables: NonEmptyList[ScatterVariableNode],
                                   scatterMethod: Option[ScatterMethod]): Checked[Set[GraphNode]] = {

    val scatterProcessingFunctionCheck = (scatterVariables.size, scatterMethod) match {
      // If we scatter over one variable only, the default processing method can handle it
      case (1, _) => (ScatterNode.DefaultScatterProcessingFunction, ScatterNode.DefaultScatterCollectionFunctionBuilder).validNelCheck
      case (_, Some(method)) => (processingFunction(method), collectingFunctionBuilder(method)).validNelCheck
      case (_, None) => "When scattering over multiple variables, a scatter method needs to be defined. See http://www.commonwl.org/v1.0/Workflow.html#WorkflowStep".invalidNelCheck
    }

    val callNode = callNodeAndNewNodes.node

    // We need to generate PBGONs for every output port of the call, so that they can be linked outside the scatter graph
    val portBasedGraphOutputNodes = callNode.outputPorts.map(op => {
      PortBasedGraphOutputNode(op.identifier, op.womType, op)
    })

    for {
      scatterFunctions <- scatterProcessingFunctionCheck
      (scatterProcessingFunction, scatterCollectingFunctionBuilder) = scatterFunctions
      // Build the scatter inner graph using the callNode, the outerGraphInputNodes, the scatterVariableNodes, and the PBGONs
      innerGraph <- Graph.validateAndConstruct(Set(callNodeAndNewNodes.node) ++ callNodeAndNewNodes.usedOuterGraphInputNodes ++ scatterVariables.toList ++ portBasedGraphOutputNodes).toEither
    } yield {
      // Build the ScatterNode - only with single scatter variable for now
      val scatterNodeWithNewNodes = ScatterNode.scatterOverGraph(innerGraph, scatterVariables.toList, scatterProcessingFunction, scatterCollectingFunctionBuilder)
      knownNodes ++ scatterNodeWithNewNodes.nodes ++ stepInputFold.generatedNodes
    }
  }

  private val CrossProductScatterProcessingFunction: ScatterProcessingFunction = { nodesAndValues: List[ScatterVariableAndValue] =>
    @tailrec
    def updateIndexLengthRec(currentList: List[ScatterVariableAndValue], combinedArraySize: Int): Int = currentList match {
      case Nil => combinedArraySize
      case ScatterVariableAndValue(variableNode, arrayValue) :: tail =>
        val indexLength = if (combinedArraySize == 0) 1 else combinedArraySize
        variableNode.withIndexLength(indexLength)
        updateIndexLengthRec(tail, indexLength * arrayValue.size)
    }

    // Reverse the list so that we start from the last variable and make our way back to the first one
    updateIndexLengthRec(nodesAndValues.reverse, 0).validNelCheck
  }

  private val NestedCrossProductScatterCollectingFunctionBuilder: ScatterCollectionFunctionBuilder = { nodesAndValues: List[ScatterVariableAndValue] =>
    (sortedShards: List[WomValue], valueType: WomArrayType) =>

      def buildCollectorArrayRec(currentList: List[Int], currentWomValues: List[WomValue], currentWomType: WomArrayType): WomArray = {
        currentList match {
          // When there's only one element left we can stop, otherwise we would end up wrapping the array twice
          case _ :: Nil | Nil => WomArray(currentWomType, currentWomValues)
          case head :: tail => val womArrays = currentWomValues.grouped(head).toList map { subArray =>
            WomArray(currentWomType, subArray)
          }
            buildCollectorArrayRec(tail, womArrays, WomArrayType(currentWomType))
        }
      }

      // Reverse sizes list, we're going to start from the most nested arrays and build up
      val res = buildCollectorArrayRec(nodesAndValues.map(_.arrayValue.size).reverse, sortedShards, valueType)
    res
  }

  def processingFunction(scatterMethod: ScatterMethod) = scatterMethod match {
    case ScatterMethod.DotProduct => ScatterNode.DefaultScatterProcessingFunction
    case ScatterMethod.FlatCrossProduct => CrossProductScatterProcessingFunction
    case ScatterMethod.NestedCrossProduct => CrossProductScatterProcessingFunction
  }

  def collectingFunctionBuilder(scatterMethod: ScatterMethod) = scatterMethod match {
    case ScatterMethod.DotProduct => ScatterNode.DefaultScatterCollectionFunctionBuilder
    case ScatterMethod.FlatCrossProduct => ScatterNode.DefaultScatterCollectionFunctionBuilder
    case ScatterMethod.NestedCrossProduct => NestedCrossProductScatterCollectingFunctionBuilder
  }
}
