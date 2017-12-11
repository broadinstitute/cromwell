package cwl

import cats.data.NonEmptyList
import cats.syntax.validated._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cwl.ScatterMethod.ScatterMethod
import cwl.WorkflowStep.WorkflowStepInputFold
import cwl.command.ParentName
import shapeless.{Coproduct, Poly1}
import wom.callable.Callable.InputDefinition
import wom.graph.CallNode.{CallNodeAndNewNodes, CallNodeBuilder, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNodePort.{OutputPort, ScatterGathererPort}
import wom.graph.ScatterNode.{ScatterCollectionFunctionBuilder, ScatterNodeBuilder, ScatterProcessingFunction, ScatterVariableAndValue}
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.types.{WomArrayType, WomMaybeEmptyArrayType, WomType}
import wom.values.{WomArray, WomValue}

import scala.annotation.tailrec

/**
  * Contains methods used to implement the scatter related logic for CWL workflows.
  * The general idea is as follows:
  * 
  * # WOM structure:
  * 
  * - Every workflow step input gets an ExpressionNode (true even for non scattered steps)
  * - Every workflow step input **being scattered over** gets its corresponding ScatterVariableNode (SVN)
  * - The ScatterNode will be depending on all the workflow step input being scattered over
  * 
  * In the scatter inner graph, will be found the CallNode, the generated PortBasedGraphOutputNodes (PBGON) for each outputs of the call,
  * the SVNs, as well as OuterGraphInputNodes (OGIN) for the the call inputs that do not point to a step input being scattered over.
  * In that case they will point to an OuterGraphInputNode (in the inner graph) that will itself reference the outer ExpressionNode
  * of the step input.
  * 
  * The InputDefinitionPointers will be pointing to either an SVN or an OGIN depending on
  * whether or not the corresponding step input is being scattered over.
  * 
  * e.g for a workflow with 2 inputs, and a single step scattering over the first input:
  *  ________________    ________________
  * |                |  |                |
  * | WorkflowInput1 |  | WorkflowInput2 |
  * |________________|  |________________|
  *          |                   |
  *  ________|_______    ________|_______
  * |                |  |                |
  * |   StepInput1   |  |   StepInput2   |
  * |________________|  |________________|
  *                |              |
  *  ______________|______________|______
  * |ScatterNode   |              |     |
  * |              |              |     |
  * |             SVN            OGIN   |
  * |              |              |     |
  * |              |              |     |
  * |              |-- CallNode --      |
  * |                     |             |
  * |                   PBGON           |
  * |___________________________________|
  * 
  * Note that the scatter node really only depends on StepInput1 (the only variable being scattered over).
  * 
  * # Scattering over multiple variables
  * 
  * When several variables are being scattered over, a scatter method is mandatory to specify how the elements must be combined.
  * Let's use the following 2 arrays as an example:
  * 
  * a1: ["one", "two"], a2: ["three", "four"]
  * 
  * * DotProduct:
  * DotProduct only works if all the arrays have the same length.
  * 
  * Shard 0: "one" - "three"
  * Shard 1: "two" - "four"
  * 
  * * CrossProduct (nested or flat):
  * Both cross product methods generate the same shards, the difference is in how they are collected at the end
  *
  * Shard 0: "one" - "three"
  * Shard 1: "one" - "four"
  * Shard 2: "two" - "three"
  * Shard 3: "two" - "four"
  * 
  * To support this, each SVN will have a method which given a shard index will return the index at which the value should be looked up in the array.
  * 
  * For dot product, this function is the same for all SVN and will be the identity[Int] function.
  * For cross product, this function depends on the number of elements in each array, which is known at runtime.
  *   When we do have this information we update each SVN with it in the ScatterProcessingFunction and returns the number of shards to be generated.
  */
object ScatterLogic {

  object ScatterVariablesPoly extends Poly1 {
    implicit def fromString: Case.Aux[String, List[String]] = at[String] { s: String => List(s) }
    implicit def fromStringList: Case.Aux[Array[String], List[String]] = at[Array[String]] { l: Array[String] => l.toList }
  }

  // Validate that the scatter expression is an ArrayLike and return the member type
  private [cwl] def scatterExpressionItemType(expressionNode: ExpressionNode) = {
    expressionNode.womType match {
      case WomArrayType(itemType) => itemType.validNelCheck // Covers maps because this is a custom unapply (see WomArrayType)
      case other => s"Cannot scatter over a non-traversable type ${other.toDisplayString}".invalidNelCheck
    }
  }

  // Build a list (potentially empty) of scatter variable nodes. Each node represents an input variable being scattered over
  def buildScatterVariables(scatter: ScatterVariables, stepInputFold: WorkflowStepInputFold, stepId: String)(implicit parentName: ParentName): Checked[List[ScatterVariableNode]] = {
    import cats.implicits._

    def buildScatterVariable(scatterVariableName: String): Checked[ScatterVariableNode] = {
      // Assume the variable is a step input (is that always true ??). Find the corresponding expression node
      val parsedScatterVariable = FileStepAndId(scatterVariableName)

      stepInputFold.stepInputMapping.get(parsedScatterVariable.id) match {
        case Some(expressionNode) =>
          // find the item type - Can't map over Either in 2.11 so convert back and forth to Validated...
          (scatterExpressionItemType(expressionNode).toValidated map { itemType =>
            // create a scatter variable node for other scattered nodes to point to
            ScatterVariableNode(WomIdentifier(scatterVariableName), expressionNode, itemType)
          }).toEither
        case None => s"Could not find a variable ${parsedScatterVariable.id} in the workflow step input to scatter over. Please make sure ${parsedScatterVariable.id} is an input of step $stepId".invalidNelCheck
      }
    }
    // Take the scatter field defining the (list of) input(s) to be scattered over
    scatter
      // Fold them so we have a List[String] no matter what (a single scattered input becomes a single element list)
      .map(_.fold(ScatterVariablesPoly))
      // If there's no scatter make it an empty list
      .getOrElse(List.empty)
      // Traverse the list to create ScatterVariableNodes that will be used later on to create the ScatterNode
      .traverse[Checked, ScatterVariableNode](buildScatterVariable)
  }

  // Build the InputDefinitionFold for an input definition for scattered steps
  def buildScatteredInputFold(inputDefinition: InputDefinition,
                              scatterVariables: List[ScatterVariableNode],
                              expressionNode: ExpressionNode,
                              callNodeBuilder: CallNodeBuilder): ErrorOr[InputDefinitionFold] = {
    // mapping = Where does this input definition get its value from (in this method it'll always be an output port, but wrapped in a Coproduct[InputDefinitionPointer])
    // outputPort = output port from which the input definition will get its values
    // newNode = potentially newly created OGIN to be added to the fold
    val (mapping, outputPort, newNode) = {
      scatterVariables.find(_.scatterExpressionNode == expressionNode) match {
        // This input is being scattered over
        case Some(scatterVariableNode) =>
          // Point the input definition to the scatter variable
          val outputPort = scatterVariableNode.singleOutputPort
          val mapping = List(inputDefinition -> Coproduct[InputDefinitionPointer](outputPort: OutputPort))
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

    val scatterGatherPortTypeFunction: WomType => WomArrayType = scatterMethod match {
      case Some(ScatterMethod.NestedCrossProduct) =>
        innerType: WomType =>
          scatterVariables.tail.foldLeft(WomArrayType(innerType))({ case (t, _) => WomArrayType(t) })
      case _ => innerType: WomType => WomArrayType(innerType)
    }

    val callNode = callNodeAndNewNodes.node

    // We need to generate PBGONs for every output port of the call, so that they can be linked outside the scatter graph
    val portBasedGraphOutputNodes = callNode.outputPorts.map(op => PortBasedGraphOutputNode(op.identifier, op.womType, op))
    
    def buildScatterNode(innerGraph: Graph, scatterProcessingFunction: ScatterProcessingFunction, scatterCollectingFunctionBuilder: ScatterCollectionFunctionBuilder) = {
      val scatterNodeBuilder = new ScatterNodeBuilder
      val outputPorts: Set[ScatterGathererPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
        scatterNodeBuilder.makeOutputPort(scatterGatherPortTypeFunction(gon.womType), gon)
      }

      val scatterNodeWithNewNodes = scatterNodeBuilder.build(innerGraph, outputPorts, scatterVariables.toList, scatterProcessingFunction, scatterCollectingFunctionBuilder)
      knownNodes ++ scatterNodeWithNewNodes.nodes ++ stepInputFold.generatedNodes
    }

    // TODO POST 2.11 Could be forcomped in 2.12 but not in 2.11 because of the Checked
    scatterProcessingFunctionCheck match {
      case Right((scatterProcessingFunction, scatterCollectingFunctionBuilder)) =>
        Graph.validateAndConstruct(Set(
          callNodeAndNewNodes.node) ++ 
          callNodeAndNewNodes.usedOuterGraphInputNodes ++ 
          scatterVariables.toList ++ 
          portBasedGraphOutputNodes
        ).map(
          buildScatterNode(_, scatterProcessingFunction, scatterCollectingFunctionBuilder)
        ).toEither
      case Left(errors) => Left(errors)
    }
  }

  /*
   * Processing function for cross products. Logic is as follows:
   * We want 1) know how many shards should be created (which this function will return)
   *         2) update each SVN with their relative index length (see ScatterVariableNode for more detail)
   *         
   * We are going to multiply all the array sizes (which will give us 1) at the end)
   * Along the way we will update the relative index length for each SVN.
   * To do so, we start at the end of the list and recurse.
   */
  private [cwl] val CrossProductScatterProcessingFunction: ScatterProcessingFunction = { nodesAndValues: List[ScatterVariableAndValue] =>
    @tailrec
    def updateIndexLengthRec(currentList: List[ScatterVariableAndValue], combinedArraySize: Int): Int = currentList match {
      case Nil => combinedArraySize
      case ScatterVariableAndValue(variableNode, arrayValue) :: tail =>
        variableNode.withRelativeIndexLength(combinedArraySize)
        updateIndexLengthRec(tail, combinedArraySize * arrayValue.size)
    }

    // Reverse the list so that we start from the last variable and make our way back to the first one
    updateIndexLengthRec(nodesAndValues.reverse, 1).validNelCheck
  }

  private [cwl] val DotProductScatterProcessingFunction: ScatterProcessingFunction = ScatterNode.DefaultScatterProcessingFunction

  // Recursively nest shard results appropriately based on the size of the scatter variable arrays
  private [cwl] val NestedCrossProductScatterCollectionFunctionBuilder: ScatterCollectionFunctionBuilder = { arraySizes: List[Int] =>
    (sortedShards: List[WomValue], valueType: WomArrayType) =>

      @tailrec
      def buildCollectorArrayRec(currentList: List[Int], currentWomValues: List[WomValue], currentWomType: WomType): WomArray = {
        /*  We stop right before the last element because it will necessarily be the number of elements in currentWomValues (provided the list is non empty).
          * e.g if currentList is (2, 3, 2), then currentWomValues will have 2 * 3 * 2 = 12 elements (this is a cross product)
          * As we recurse, we keep grouping currentWomValues by the numbers in currentList:
          * 
          * group by 2: [[_, _], [_, _], [_, _], [_, _], [_, _], [_, _]]
          * group by 3: [
          *   [[_, _], [_, _], [_, _]],
          *   [[_, _], [_, _], [_, _]]
          * ]
          * This list has necessarily 2 elements already because 12 / 2 / 3 = 2. So we can make the final WomArray out of that
        */
        currentList match {
          case _ :: Nil | Nil => 
            WomArray(WomMaybeEmptyArrayType(currentWomType), currentWomValues)
          case head :: tail => 
            val arrayType = WomMaybeEmptyArrayType(currentWomType)
            val womArrays = currentWomValues.grouped(head).toList map { WomArray(arrayType, _) }
            buildCollectorArrayRec(tail, womArrays, arrayType)
        }
      }
      
      def mostInnerType(tpe: WomType): WomType = tpe match {
        case WomArrayType(t) => mostInnerType(t)
        case other => other
      }
      
      def buildEmptyArrays(arraySizeList: List[Int], womType: WomType, womValues: List[WomValue]): WomArray = arraySizeList match {
        case Nil => womType match {
          case arrayType: WomArrayType => WomArray(arrayType, womValues)
          case _ => 
            // This would mean arraySizeList was empty to begin with (otherwise we would have wrapped womType in a WomArrayType at least once)
            // which should be impossible since it would mean there was no scatter variables, so we shouldn't even be here
            throw new RuntimeException("Programmer error ! We should not be collecting scatter nodes if there was no scatter !")
        }
        case 0 :: tail =>
          buildEmptyArrays(tail, WomArrayType(womType), List.empty)
        case head :: tail =>
          val arrayType = womType match {
            case array: WomArrayType => array
            case nonArray => WomArrayType(nonArray)
          }
          val womArrays = (0 until head).toList map { _ => WomArray(arrayType, womValues) }
          buildEmptyArrays(tail, WomArrayType(womType), womArrays)
      }
      
      /*
       * There are 2 distinct cases.
       * If we have no shards, it means that at least one of the scatter expressions evaluated to an empty array.
       * In that case we want to create an array structure of empty arrays reflecting the scatter array sizes.
       * e.g:
       * 
       * scatter array sizes   |   collected result
       *        [0]                       []
       *       [0, 2]                     [] 
       *       [2, 0]                   [[], []]
       *      [2, 0, 3]                 [[], []]
       *      [2, 3, 0]                [[[], [], []], [[], [], []]] 
       *
       * Note that as soon as we reach a size 0, we return an empty array.
       * 
       * If we have shards, we successively group the shards list to build arrays of the right size.
       * 
       * In both cases we reverse the list so that we can left recurse while still building the result starting from the most nested arrays
       */
      if (sortedShards.isEmpty) {
        buildEmptyArrays(arraySizes.reverse, mostInnerType(valueType), List.empty)
      } else {
        // 
        buildCollectorArrayRec(arraySizes.reverse, sortedShards, mostInnerType(valueType))
      }
  }

  // select a processing function based on the scatter method
  private def processingFunction(scatterMethod: ScatterMethod) = scatterMethod match {
    case ScatterMethod.DotProduct => DotProductScatterProcessingFunction
      // Both cross product methods use the same processing function, the difference is in how we collect results (see collectingFunctionBuilder)
    case ScatterMethod.FlatCrossProduct => CrossProductScatterProcessingFunction
    case ScatterMethod.NestedCrossProduct => CrossProductScatterProcessingFunction
  }

  // select a collecting function builder based on the scatter method
  private def collectingFunctionBuilder(scatterMethod: ScatterMethod) = scatterMethod match {
      // dot product and flat cross product output a flat array, which is the default behavior
    case ScatterMethod.DotProduct => ScatterNode.DefaultScatterCollectionFunctionBuilder
    case ScatterMethod.FlatCrossProduct => ScatterNode.DefaultScatterCollectionFunctionBuilder
      // nested cross product uses a special collecting function to build nested arrays
    case ScatterMethod.NestedCrossProduct => NestedCrossProductScatterCollectionFunctionBuilder
  }
}
