package wdl

import cats.data.Validated.Valid
import cats.syntax.validated._
import common.collections.EnhancedCollections._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.AstTools.{EnhancedAstNode, VariableReference}
import wdl.Declaration.{GraphOutputDeclarationNode, InputDeclarationNode, IntermediateValueDeclarationNode}
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExposedExpressionNode


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

    // We want:
    // - Nodes that this node references
    // - The immediate ancestor
    // - Anything that our children sre downstream from
    // But because our children's upstream might also include these (which we don't want), filter out:
    // - This
    // - Any other WdlGraphNode descendants of this
    (referencedNodes ++ closestScopedAncestor.toSeq ++ childGraphNodes.flatMap(_.upstream)).toSet - this -- descendants.filterByType[WdlGraphNode]
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
  private[wdl] def buildWomGraph(from: Scope, includeGraphNodes: Set[GraphNode], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[Graph] = {

    final case class FoldState(nodes: Set[GraphNode], availableInputs: Map[String, GraphNodePort.OutputPort])

    val initialFoldState = FoldState(
      Set.empty,
      includeGraphNodes.flatMap { gn =>
        gn.outputPorts map { p => p.name -> p }
      }.toMap
    )

    def foldFunction(acc: ErrorOr[FoldState], node: WdlGraphNode): ErrorOr[FoldState] = acc flatMap { goodAcc => buildNode(goodAcc, node) }

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
      case wdlCall: WdlCall => WdlCall.buildWomNodeAndInputs(wdlCall, acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { cnani: CallNodeAndNewNodes =>
        foldInGeneratedNodeAndNewInputs(acc, cnani.node.localName + ".")(cnani)
      }

      case decl: DeclarationInterface => Declaration.buildWdlDeclarationNode(decl, acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { wdlDeclNode =>
        // As with GeneratedNodeAndNewInputs, we might have made some new OuterGraphInputNodes to build this DeclarationNode, so
        // make sure they get included:
        val declNode = wdlDeclNode.toGraphNode
        val newOgins: Set[OuterGraphInputNode] = declNode.upstreamOuterGraphInputNodes
        val newOginOutputs = newOgins.map(_.nameToPortMapping)

        // We only want to include the output port as an input to other GraphNodes in the FoldState if it's not a Workflow Output:
        val availableOutputPort: Option[(String, OutputPort)] = wdlDeclNode match {
          case InputDeclarationNode(graphInputNode) => Option(declNode.localName -> graphInputNode.singleOutputPort)
          case IntermediateValueDeclarationNode(expressionNode) => Option(declNode.localName -> expressionNode.singleExpressionOutputPort)
          case GraphOutputDeclarationNode(_) => None
        }

        FoldState(acc.nodes + declNode ++ newOgins, acc.availableInputs ++ newOginOutputs ++ availableOutputPort)
      }

      case scatter: Scatter =>
        Scatter.womScatterNode(scatter, acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { foldInGeneratedNodeAndNewInputs(acc, "")(_) }
      case ifBlock: If =>
        If.womConditionalNode(ifBlock, acc.availableInputs, outerLookup, preserveIndexForOuterLookups) map { foldInGeneratedNodeAndNewInputs(acc, "")(_) }

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
