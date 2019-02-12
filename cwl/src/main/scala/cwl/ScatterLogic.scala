package cwl

import cats.data.NonEmptyList
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import cwl.ScatterMethod.ScatterMethod
import shapeless.Poly1
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNodePort.ScatterGathererPort
import wom.graph.ScatterNode._
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.types.{WomArrayType, WomMaybeEmptyArrayType, WomType}
import wom.values.{WomArray, WomValue}

import scala.annotation.tailrec

/**
  * Contains methods used to implement the scatter related logic for CWL workflows.
  * The general idea is as follows:
  *
  * See WorkflowStep for an example of WOM graph for scatters.
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
  case class ScatterFunctions(processing: ScatterProcessingFunction, collection: ScatterCollectionFunctionBuilder)

  object ScatterVariablesPoly extends Poly1 {
    implicit def fromString: Case.Aux[String, List[String]] = at[String] { s: String => List(s) }
    implicit def fromStringList: Case.Aux[Array[String], List[String]] = at[Array[String]] { l: Array[String] => l.toList }
  }

  // Validate that the scatter expression is an ArrayLike and return the member type
  private [cwl] def scatterExpressionItemType(expressionNode: ExpressionNode) = {
    expressionNode.womType match {
      case WomArrayType(itemType) => itemType.validNelCheck // Covers maps because this is a custom unapply (see WomArrayType)
      case other => s"Cannot scatter over a non-traversable type ${other.stableName}".invalidNelCheck
    }
  }

  def scatterGatherPortTypeFunction(scatterMethod: Option[ScatterMethod], scatterVariables: NonEmptyList[_]): WomType => WomArrayType = scatterMethod match {
    case Some(ScatterMethod.NestedCrossProduct) =>
      innerType: WomType =>
        scatterVariables.tail.foldLeft(WomArrayType(innerType))({ case (t, _) => WomArrayType(t) })
    case _ => innerType: WomType => WomArrayType(innerType)
  }

  // Build a map (potentially empty) of scatterered input steps with their corresponding SVN node.
  def buildScatterVariableNodes(scatter: ScatterVariables, stepInputMappings: Map[WorkflowStepInput, ExpressionNode], stepId: String): Checked[Map[WorkflowStepInput, ScatterVariableNode]] = {
    import cats.implicits._

    def buildScatterVariable(scatterVariableName: String): ErrorOr[(WorkflowStepInput, ScatterVariableNode)] = {
      // Assume the variable is a step input (is that always true ??). Find the corresponding expression node

      stepInputMappings.find({ case (stepInput, _) => stepInput.id == scatterVariableName }) match {
        case Some((stepInput, expressionNode)) =>
          scatterExpressionItemType(expressionNode).toValidated map { itemType =>
            // create a scatter variable node for other scattered nodes to point to
            stepInput -> ScatterVariableNode(WomIdentifier(scatterVariableName), expressionNode, itemType)
          }
        case None => s"Could not find a variable $scatterVariableName in the workflow step input to scatter over. Please make sure $scatterVariableName is an input of step $stepId".invalidNel
      }
    }
    // Take the scatter field defining the (list of) input(s) to be scattered over
    scatter
      // Fold them so we have a List[String] no matter what (a single scattered input becomes a single element list)
      .map(_.fold(ScatterVariablesPoly))
      // If there's no scatter make it an empty list
      .getOrElse(List.empty)
      // Traverse the list to create ScatterVariableNodes that will be used later on to create the ScatterNode
      .traverse(buildScatterVariable)
      .toEither
      .map(_.toMap)
  }

  // Prepare the nodes to be returned if this call is being scattered
  def buildScatterNode(callNodeAndNewNodes: CallNodeAndNewNodes,
                       scatterVariableNodes: NonEmptyList[ScatterVariableNode],
                       ogins: Set[OuterGraphInputNode],
                       stepExpressionNodes: Set[ExpressionNode],
                       scatterMethod: Option[ScatterMethod]): Checked[ScatterNode] = {

    val scatterProcessingFunctionCheck = (scatterVariableNodes.size, scatterMethod) match {
      // If we scatter over one variable only, the default processing method can handle it
      case (1, _) => ScatterFunctions(ScatterNode.DefaultScatterProcessingFunction, ScatterNode.DefaultScatterCollectionFunctionBuilder).validNelCheck
      case (_, Some(method)) => ScatterFunctions(processingFunction(method), collectingFunctionBuilder(method)).validNelCheck
      case (_, None) => "When scattering over multiple variables, a scatter method needs to be defined. See http://www.commonwl.org/v1.0/Workflow.html#WorkflowStep".invalidNelCheck
    }

    val callNode = callNodeAndNewNodes.node

    // We need to generate PBGONs for every output port of the call, so that they can be linked outside the scatter graph
    val portBasedGraphOutputNodes = callNode.outputPorts.map(op => PortBasedGraphOutputNode(op.identifier, op.womType, op))

    def buildScatterNode(innerGraph: Graph, scatterFunctions: ScatterFunctions) = {
      val scatterNodeBuilder = new ScatterNodeBuilder
      val outputPorts: Set[ScatterGathererPort] = innerGraph.nodes.collect { case gon: PortBasedGraphOutputNode =>
        scatterNodeBuilder.makeOutputPort(scatterGatherPortTypeFunction(scatterMethod, scatterVariableNodes)(gon.womType), gon)
      }

      scatterNodeBuilder.build(innerGraph, outputPorts, scatterVariableNodes.toList, scatterFunctions.processing, scatterFunctions.collection)
    }
    
    for {
      scatterProcessingFunctionAndBuilder <- scatterProcessingFunctionCheck
      graph <- Graph.validateAndConstruct(
        Set(callNodeAndNewNodes.node) ++
          ogins ++
          stepExpressionNodes ++
          scatterVariableNodes.toList ++
          portBasedGraphOutputNodes
      ).toEither
      scatterNode = buildScatterNode(graph, scatterProcessingFunctionAndBuilder)
    } yield scatterNode.node
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
