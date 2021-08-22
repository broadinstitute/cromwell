package wom.views

import java.util.concurrent.atomic.AtomicInteger

import cats.syntax.all._
import cats.instances.list._
import cats.Monoid
import wom.callable.ExecutableCallable
import wom.graph.GraphNodePort.{ConditionalOutputPort, OutputPort, ScatterGathererPort}
import wom.graph.expression.ExpressionNode
import wom.graph._
import wom.types.WomType
import wom.views.GraphPrint._

final class GraphPrint(executableCallable: ExecutableCallable) {

  def dotString: String =
    WorkflowDigraph(
      workflowName = executableCallable.name,
      digraph = listAllGraphNodes(executableCallable.graph, new AtomicInteger(0), Map.empty),
    ).dotString

  // A "monoid" is just a fancy way of saying "thing you can add together".
  // The cats library makes it easy to turn case classes into Monoids - ie into "things you can add together".
  // Adding one NodesAndLinks collection to another is the same as adding together the nodes set, and the links set, of each,
  //  and making a new NodesAndLinks out of the combined results.
  // This line uses cats to allow NodesAndLinks to be added together, which is then used later on by calls to "foldMap".
  implicit val nodeAndLinkMonoid: Monoid[NodesAndLinks] = cats.derived.MkMonoid[NodesAndLinks]

  /*
   * NB: We don't care that the counter is "atomic", this is just a way to have a "reference-to-a-mutable-Int" to pass around
   */
  private def listAllGraphNodes(graph: Graph,
                                clusterCounter: AtomicInteger,
                                availableScatterVariables: Map[ScatterVariableNode, DotScatterVariableNode]): NodesAndLinks = {

    graph.nodes.toList.filter(worthDisplaying).foldMap {
      case ccn: CommandCallNode => NodesAndLinks(Set(DotCallNode(ccn)), upstreamLinks(ccn, DotCallNode(ccn), availableScatterVariables))
      case scn: WorkflowCallNode => NodesAndLinks(Set(DotSubworkflowCallNode(scn)), upstreamLinks(scn, DotSubworkflowCallNode(scn), availableScatterVariables))
      case s: ScatterNode => handleScatter(s, clusterCounter, availableScatterVariables)
      case c: ConditionalNode => handleConditional(c, clusterCounter, availableScatterVariables)
      case _ => nodeAndLinkMonoid.empty
    }
  }

  def handleScatter(scatterNode: ScatterNode,
                    clusterCounter: AtomicInteger,
                    knownScatterVariables: Map[ScatterVariableNode, DotScatterVariableNode]): NodesAndLinks = {
    val clusterNumber = clusterCounter.getAndIncrement()
    val id = s"cluster_$clusterNumber"

    def handleScatterVariableNode(node: ScatterVariableNode): NodesAndLinks = {
      val dotNode = DotScatterVariableNode.apply(node, clusterNumber)
      val links = upstreamLinks(node.linkToOuterGraph.graphNode, dotNode, Map.empty)
      NodesAndLinks(Set(dotNode), links)
    }

    val scatterExpressionNodesAndLinks = scatterNode.scatterVariableNodes.foldMap(handleScatterVariableNode)

    val madeScatterVariableNodes = scatterNode.scatterVariableNodes.map { svn =>
      val link = scatterExpressionNodesAndLinks.nodes.collectFirst {
        case dsvn: DotScatterVariableNode if dsvn.valueName == svn.localName => dsvn
      }

      svn -> link.get
    }.toMap

    val innerGraphNodesAndLinks = listAllGraphNodes(scatterNode.innerGraph, clusterCounter, knownScatterVariables ++ madeScatterVariableNodes)

    NodesAndLinks(
      nodes = Set(DotScatterNode(id, innerGraphNodesAndLinks.nodes ++ scatterExpressionNodesAndLinks.nodes)),
      links = scatterExpressionNodesAndLinks.links ++ innerGraphNodesAndLinks.links
    )
  }

  def handleConditional(conditionalNode: ConditionalNode,
                        clusterCounter: AtomicInteger,
                        knownScatterVariables: Map[ScatterVariableNode, DotScatterVariableNode]): NodesAndLinks = {
    val clusterNumber = clusterCounter.getAndIncrement()
    val id = s"cluster_$clusterNumber"

    def handleConditionalExpressionNode(conditionalExpression: ExpressionNode): NodesAndLinks = {
      val dotNode = DotConditionalExpressionNode(conditionalExpression, clusterNumber)
      val links = upstreamLinks(conditionalExpression, dotNode, Map.empty)
      NodesAndLinks(Set(dotNode), links)
    }

    val conditionalExpressionNodesAndLinks = handleConditionalExpressionNode(conditionalNode.conditionExpression)

    val innerGraphNodesAndLinks = listAllGraphNodes(conditionalNode.innerGraph, clusterCounter, knownScatterVariables)

    NodesAndLinks(
      nodes = Set(DotConditionalNode(id, innerGraphNodesAndLinks.nodes ++ conditionalExpressionNodesAndLinks.nodes)) ,
      links = innerGraphNodesAndLinks.links ++ conditionalExpressionNodesAndLinks.links
    )
  }

}

object GraphPrint {
  final case class DotLink(from: DotNode, to: DotNode) {
    def dotString = s"${from.id} -> ${to.id}"
  }

  final case class NodesAndLinks(nodes: Set[DotNode], links: Set[DotLink])

  final case class WorkflowDigraph(workflowName: String, digraph: NodesAndLinks) {
    def dotString: String =
      s"""|digraph $workflowName {
          |  #rankdir=LR;
          |  compound=true;
          |
          |  # Links
          |  ${digraph.links.toList.flatMap(_.dotString.linesIterator).mkString(System.lineSeparator + "  ")}
          |
          |  # Nodes
          |  ${digraph.nodes.toList.flatMap(_.dotString.linesIterator).mkString(System.lineSeparator + "  ")}
          |}""".stripMargin
  }

  sealed trait DotNode { def id: String; def dotString: String }

  final case class DotCallNode(callName: String) extends DotNode {
    override def id: String = s"CALL_$callName"
    def dotString = s"""$id [label="call $callName"]"""
  }
  object DotCallNode {
    def apply(ccn: CommandCallNode): DotCallNode = DotCallNode(ccn.localName)
  }

  final case class DotSubworkflowCallNode(callName: String) extends DotNode {
    override def id: String = s"CALL_$callName"
    def dotString = s"""$id [label="call $callName";shape="oval";peripheries=2]"""
  }
  object DotSubworkflowCallNode {
    def apply(scn: WorkflowCallNode): DotSubworkflowCallNode = DotSubworkflowCallNode(scn.localName)
  }

  final case class DotScatterVariableNode(womType: WomType, valueName: String, clusterNumber: Int) extends DotNode {
    override def id: String = s"SCATTER_${clusterNumber}_VARIABLE_$valueName"
    def dotString = s"""$id [shape="hexagon" label="scatter over ${womType.friendlyName} as $valueName"]"""
  }
  object DotScatterVariableNode {
    def apply(svn: ScatterVariableNode, clusterNumber: Int): DotScatterVariableNode = DotScatterVariableNode(svn.womType, svn.identifier.localName.value, clusterNumber)
  }

  final case class DotConditionalExpressionNode(womType: WomType, expressionString: String, clusterNumber: Int) extends DotNode {
    override def id: String = s"CONDITIONAL_${clusterNumber}_EXPRESSION"
    def dotString = s"""$id [shape="hexagon" label="if ($expressionString)" style="dashed" ]"""
  }
  object DotConditionalExpressionNode {
    def apply(en: ExpressionNode, clusterNumber: Int): DotConditionalExpressionNode = DotConditionalExpressionNode(en.womType, escapeQuotes(en.womExpression.sourceString), clusterNumber)
  }

  final case class DotScatterNode(id: String, nodes: Set[DotNode]) extends DotNode {
    override def dotString: String =
      s"""subgraph $id {
          |  style="filled,solid";
          |  fillcolor=white;
          |  ${nodes.toList.flatMap(_.dotString.linesIterator).mkString(System.lineSeparator() + "  ")}
          |}""".stripMargin
  }

  final case class DotConditionalNode(id: String, nodes: Set[DotNode]) extends DotNode {
    override def dotString: String =
      s"""subgraph $id {
         |  style="filled,dashed";
         |  fillcolor=white;
         |  ${nodes.toList.flatMap(_.dotString.linesIterator).mkString(System.lineSeparator() + "  ")}
         |}""".stripMargin
  }

  def upstreamLinks(originNode: GraphNode,
                    origin: DotNode,
                    availableScatterVariables: Map[ScatterVariableNode, DotScatterVariableNode]): Set[DotLink] = {
    def relevantAsUpstream(nodeToLink: GraphNode): Set[DotNode] = nodeToLink match {
      case ccn: CommandCallNode => Set(DotCallNode(ccn))
      case scn: WorkflowCallNode => Set(DotSubworkflowCallNode(scn))
      case en: ExpressionNode => upstreamLinksforNode(en)
      case svn: ScatterVariableNode => Set(availableScatterVariables(svn))
      case ogin: OuterGraphInputNode => upstreamPortToRelevantNodes(ogin.linkToOuterGraph)

      case _ => Set.empty[DotNode]
    }

    def upstreamPortToRelevantNodes(p: OutputPort) = p match {
      case gatherPort: ScatterGathererPort => relevantAsUpstream(gatherPort.outputToGather.singleUpstreamNode)
      case conditionalOutputPort: ConditionalOutputPort => relevantAsUpstream(conditionalOutputPort.outputToExpose.singleUpstreamNode)
      case other =>
        relevantAsUpstream(other.graphNode)
    }

    def upstreamLinksforNode(n: GraphNode) = n.upstreamPorts flatMap upstreamPortToRelevantNodes
    upstreamLinksforNode(originNode).map(DotLink(_, origin))
  }

  def hasCallAncestor(g: GraphNode): Boolean = g.upstreamAncestry.exists(_.isInstanceOf[CommandCallNode])

  def escapeQuotes(s: String): String = s.replace("\"", "\\\"")

  def worthDisplaying(node: GraphNode): Boolean = node match {
    case _: CommandCallNode => true
    case _: WorkflowCallNode => true
    case s: ScatterNode => s.innerGraph.nodes.exists(worthDisplaying)
    case c: ConditionalNode => c.innerGraph.nodes.exists(worthDisplaying)
    case _ => false
  }
}
