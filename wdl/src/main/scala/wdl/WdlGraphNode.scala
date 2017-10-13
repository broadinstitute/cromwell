package wdl

import cats.data.Validated.Valid
import cats.syntax.validated._
import lenthall.collections.EnhancedCollections._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
import wdl.AstTools.{EnhancedAstNode, VariableReference}
import wom.graph._
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort.OutputPort


sealed trait WdlGraphNode extends Scope {

  /**
    * The set of all graph nodes which are (transitively) upstream from this one.
    */
  final lazy val upstreamAncestry: Set[WdlGraphNode] = WdlGraphNode.calculateUpstreamAncestry(Set.empty, this)

  def referencedNodes: Iterable[WdlGraphNode]

  /**
    * The set of all graph nodes which are a single step upstream from this one.
    */
  final lazy val upstream: Set[WdlGraphNode] = {
    // If we are inside the scope of another graph node (i.e. via document element ancestry), then
    // that is also upstream of us.
    val closestScopedAncestor = ancestry.collectFirst({
      case ancestor: WdlGraphNode => ancestor
    })

    (referencedNodes ++ closestScopedAncestor.toSeq).toSet
  }

  final lazy val downstream: Set[WdlGraphNode] = {
    for {
      node <- namespace.descendants.collect({
        case n: WdlGraphNode if n.fullyQualifiedName != fullyQualifiedName => n
      })
      if node.upstream.contains(this)
    } yield node
  }

  def isUpstreamFrom(other: WdlGraphNode): Boolean = {
    other.upstreamAncestry.contains(this) || (other.childGraphNodes exists isUpstreamFrom)
  }
}

object WdlGraphNode {
  private def calculateUpstreamAncestry(currentSet: Set[WdlGraphNode], graphNode: WdlGraphNode): Set[WdlGraphNode] = {
    val setWithUpstream = currentSet ++ graphNode.upstream
    val updatesNeeded = graphNode.upstream -- currentSet
    updatesNeeded.foldLeft(setWithUpstream)(calculateUpstreamAncestry)
  }

  /**
    * Build a Wom Graph from a Scope in WDL.
    * - Looks at the children of the Scope and converts them into GraphNodes.
    * - Adds in any extras from 'includeGraphNode' which the call knows should also be in the Graph (eg the Scatter variable's InputGraphNode which this method doesn't know exists)
    * - Builds them all together into a Graph
    */
  private[wdl] def buildWomGraph(from: Scope, includeGraphNodes: Set[GraphNode], outerLookup: Map[String, GraphNodePort.OutputPort]): ErrorOr[Graph] = {

    final case class FoldState(nodes: Set[GraphNode], availableInputs: Map[String, GraphNodePort.OutputPort])

    val initialFoldState = FoldState(
      Set.empty,
      includeGraphNodes.flatMap { gn =>
        gn.outputPorts map { p => p.name -> p }
      }.toMap
    )

    def foldFunction(acc: ErrorOr[FoldState], node: WdlGraphNode): ErrorOr[FoldState] = acc flatMap { goodAcc =>  buildNode(goodAcc, node) }

    def foldInGeneratedNodeAndNewInputs(acc: FoldState, outputPortPrefix: String)(gnani: GeneratedNodeAndNewNodes): FoldState = {
      // The output ports from the new node
      val newCallOutputPorts = (gnani.node.outputPorts map { p => s"$outputPortPrefix${p.name}" -> p }).toMap
      // The output ports from newly created GraphInputNodes:
      val newInputOutputPorts = (gnani.newInputs map { i => i.localName -> i.singleOutputPort }).toMap

      FoldState(acc.nodes + gnani.node ++ gnani.newInputs ++ gnani.newExpressions, acc.availableInputs ++ newCallOutputPorts ++ newInputOutputPorts)
    }

    def buildNode(acc: FoldState, node: WdlGraphNode): ErrorOr[FoldState] = node match {
      case wdlCall: WdlCall => WdlCall.buildWomNodeAndInputs(wdlCall, acc.availableInputs, outerLookup) map { case cnani @ CallNodeAndNewNodes(call, _, _) =>
        foldInGeneratedNodeAndNewInputs(acc, call.localName + ".")(cnani)
      }

      case decl: DeclarationInterface => Declaration.buildWomNode(decl, acc.availableInputs, outerLookup) map { declNode =>
        FoldState(acc.nodes + declNode.toGraphNode, acc.availableInputs ++ declNode.singleOutputPort.collect { case sop: OutputPort => declNode.toGraphNode.localName -> sop })
      }

      case scatter: Scatter => Scatter.womScatterNode(scatter, acc.availableInputs) map { foldInGeneratedNodeAndNewInputs(acc, "")(_) }
      case ifBlock: If => If.womConditionalNode(ifBlock, acc.availableInputs) map { foldInGeneratedNodeAndNewInputs(acc, "")(_) }

      case _ => s"Cannot process WdlGraphNodes of type ${node.getClass.getSimpleName} yet!".invalidNel
    }

    val nodeList = from.childGraphNodesSorted
    val nodeAccumulator: ErrorOr[FoldState] = nodeList.foldLeft[ErrorOr[FoldState]](Valid(initialFoldState))(foldFunction)

    def outerLinkInputs(nodes: Set[GraphNode]): Set[OuterGraphInputNode] = nodes flatMap {
      // NB: this curious type annotation just gives intelliJ a hand:
      _.inputPorts.map(_.upstream.graphNode).filterByType[OuterGraphInputNode]: Set[OuterGraphInputNode]
    }

    def withDefaultOutputs(g: Graph): Graph = if (g.nodes.exists(_.isInstanceOf[GraphOutputNode])) { g } else {
      Graph(g.nodes.union((g.nodes collect {
          case node: CallNode => node.outputPorts.map(op => {
            val identifier = node.identifier.combine(op.name)
            PortBasedGraphOutputNode(identifier, op.womType, op)
          })
          case node: ScatterNode => node.outputMapping.map(op => {
            PortBasedGraphOutputNode(op.identifier, op.womType, op)
          })
          case node: ConditionalNode => node.conditionalOutputPorts.map(op => {
            PortBasedGraphOutputNode(op.identifier, op.womType, op)
          })
        }).flatten))
    }

    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    for {
      foldState <- nodeAccumulator
      graphNodes = foldState.nodes
      outerLinks = outerLinkInputs(graphNodes)
      g <- Graph.validateAndConstruct(graphNodes ++ outerLinks ++ includeGraphNodes)
      withOutputs = withDefaultOutputs(g)
    } yield withOutputs
  }
}

trait WdlGraphNodeWithUpstreamReferences extends WdlGraphNode {
  def upstreamReferences: Iterable[VariableReference]

  // If we have variable reference to other graph nodes, then they are upstream from us.
  override final def referencedNodes = for {
      variable <- upstreamReferences
      node <- resolveVariable(variable.terminal.sourceString)
      if node.fullyQualifiedNameWithIndexScopes != fullyQualifiedNameWithIndexScopes
    } yield node
}

trait WdlGraphNodeWithInputs extends WdlGraphNode {
  def inputMappings: Map[String, WdlExpression]

  override final def referencedNodes = for {
    expr <- inputMappings.values
    variable <- expr.variableReferences
    scope <- parent
    node <- scope.resolveVariable(variable.terminal.sourceString)
  } yield node
}
