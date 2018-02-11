package wdl.transforms.draft2.wdlom2wom

import cats.data.Validated.Valid
import cats.syntax.validated._
import common.collections.EnhancedCollections._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.draft2.model._
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.ExposedExpressionNode
import wom.graph._
import wom.transforms.WomCallNodeMaker.ops._
import wom.transforms.WomConditionalNodeMaker.ops._
import wom.transforms.WomScatterNodeMaker.ops._
import wom.transforms.WomGraphMaker

object WdlDraft2WomGraphMaker extends WomGraphMaker[Scope] {

  /**
    * Build a Wom Graph from a Scope in WDL.
    * - Looks at the children of the Scope and converts them into GraphNodes.
    * - Adds in any extras from 'includeGraphNode' which the call knows should also be in the Graph (eg the Scatter variable's InputGraphNode which this method doesn't know exists)
    * - Builds them all together into a Graph
    */
  override def toWomGraph(scope: Scope, includeGraphNodes: Set[GraphNode], outerLookup: Map[String, OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[Graph] = {

    final case class FoldState(nodes: Set[GraphNode], availableInputs: Map[String, GraphNodePort.OutputPort])

    val initialFoldState = FoldState(
      Set.empty,
      includeGraphNodes.flatMap { gn =>
        gn.outputPorts map { p => p.name -> p }
      }.toMap
    )

    def foldFunction(acc: ErrorOr[FoldState], node: WdlGraphNode): ErrorOr[FoldState] = acc flatMap { goodAcc =>
      buildNode(goodAcc, node).leftMap(errors => errors.map(s"Unable to build WOM node for ${node.getClass.getSimpleName} '${node.womIdentifier.localName.value}': " + _))
    }

    def foldInGeneratedNodeAndNewInputs(acc: FoldState, outputPortPrefix: String)(gnani: GeneratedNodeAndNewNodes): FoldState = {
      // The output ports from the new node
      val newCallOutputPorts = (gnani.node.outputPorts map { p => s"$outputPortPrefix${p.name}" -> p }).toMap
      // The output ports from newly created GraphInputNodes:
      val newInputOutputPorts = (gnani.newInputs map { i => i.localName -> i.singleOutputPort }).toMap
      val usedOuterGraphInputNodes = gnani.usedOuterGraphInputNodes
      val newOginOutputs = usedOuterGraphInputNodes.map(_.nameToPortMapping)

      // To make our new list of Nodes in this graph, we add:
      // - All the existing Nodes
      // - The newly generated Node
      // - Any new inputs made by the generated Node
      // - Any OGINs that were generated during the creation of the new Node
      // - Any new Expression Nodes that were made during the generation of the new Node.
      // To make our new list of available inputs, we add:
      // - All the existing ports in the FoldState
      // - Any new call outputs from the generated Node
      // - The outputs from any new input Nodes we had to make (so that other Nodes don't have to recreate them)
      FoldState(
        nodes = acc.nodes + gnani.node ++ gnani.newInputs ++ usedOuterGraphInputNodes ++ gnani.newExpressions,
        availableInputs = acc.availableInputs ++ newCallOutputPorts ++ newInputOutputPorts ++ newOginOutputs
      )
    }

    def buildNode(acc: FoldState, node: WdlGraphNode): ErrorOr[FoldState] = node match {

      case wdlCall: WdlCall =>
        wdlCall.toWomCallNode(acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { cnani: CallNodeAndNewNodes =>
          foldInGeneratedNodeAndNewInputs(acc, cnani.node.localName + ".")(cnani)
        }

      case decl: DeclarationInterface => Declaration.buildWdlDeclarationNode(decl, acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { wdlDeclNode =>
        // As with GeneratedNodeAndNewInputs, we might have made some new OuterGraphInputNodes to build this DeclarationNode, so
        // make sure they get included:
        val declNode = wdlDeclNode.toGraphNode
        val newOgins: Set[OuterGraphInputNode] = declNode.upstreamOuterGraphInputNodes
        val newOginOutputs = newOgins.map(_.nameToPortMapping)

        // Add the output port, but only if the value isn't already available as an input from somewhere else
        val availableOutputPort: Option[(String, OutputPort)] = if(!acc.availableInputs.contains(wdlDeclNode.localName)) {
          Option(wdlDeclNode.localName -> wdlDeclNode.singleOutputPort)
        } else {
          None
        }

        FoldState(acc.nodes + declNode ++ newOgins, acc.availableInputs ++ newOginOutputs ++ availableOutputPort)
      }

      case scatter: Scatter =>
        scatter.toWomScatterNode(acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { foldInGeneratedNodeAndNewInputs(acc, "")(_) }
      case ifBlock: If =>
        ifBlock.toWomConditionalNode(acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { foldInGeneratedNodeAndNewInputs(acc, "")(_) }

      case _ => s"Cannot process WdlGraphNodes of type ${node.getClass.getSimpleName} yet!".invalidNel
    }

    val nodeList = scope.childGraphNodesSorted
    val nodeAccumulator: ErrorOr[FoldState] = nodeList.foldLeft[ErrorOr[FoldState]](Valid(initialFoldState))(foldFunction)

    def outerLinkInputs(nodes: Set[GraphNode]): Set[OuterGraphInputNode] = nodes flatMap {
      // NB: this curious type annotation just gives intelliJ a hand:
      _.inputPorts.map(_.upstream.graphNode).filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode]
    }

    // Default outputs are added if we're:
    // - A scatter block
    // - An if block
    // - (NB: top level workflows are already given wildcard outputs in the WDL Workflow building phase)
    def withDefaultOutputs(g: Graph): Graph = scope match {
      case _: If | _: Scatter =>
        Graph(g.nodes.union((g.nodes collect {
          case node: CallNode => node.outputPorts.map(op => {
            val identifier = node.identifier.combine(op.name)
            PortBasedGraphOutputNode(identifier, op.womType, op)
          })
          case node: ExposedExpressionNode => node.outputPorts.map(op => {
            PortBasedGraphOutputNode(WomIdentifier(op.name), op.womType, op)
          })
          case node: ScatterNode => node.outputMapping.map(op => {
            PortBasedGraphOutputNode(op.identifier, op.womType, op)
          })
          case node: ConditionalNode => node.conditionalOutputPorts.map(op => {
            PortBasedGraphOutputNode(op.identifier, op.womType, op)
          })
        }).flatten))
      case _ => g
    }

    import common.validation.ErrorOr.ShortCircuitingFlatMap
    for {
      foldState <- nodeAccumulator
      graphNodes = foldState.nodes
      outerLinks = outerLinkInputs(graphNodes)
      g <- Graph.validateAndConstruct(graphNodes ++ outerLinks ++ includeGraphNodes)
      withOutputs = withDefaultOutputs(g)
    } yield withOutputs
  }
}
