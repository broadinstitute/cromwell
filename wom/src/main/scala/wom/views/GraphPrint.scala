package wom.views

import java.util.concurrent.atomic.AtomicInteger

import cats.implicits._
import cats.Monoid
import wom.callable.ExecutableCallable
import wom.graph.GraphNodePort.{OutputPort, ScatterGathererPort}
import wom.graph.expression.{ExposedExpressionNode, ExpressionNode}
import wom.graph._
import wom.types.WomType
import wom.views.GraphPrint._

final class GraphPrint(executableCallable: ExecutableCallable) {

  def dotString = WorkflowDigraph(executableCallable.name, listAllGraphNodes(executableCallable.graph, new AtomicInteger(0), Map.empty)).dotString

  implicit val nodeAndLinkMonoid: Monoid[NodesAndLinks] = cats.derived.MkMonoid[NodesAndLinks]

  /*
   * NB: We don't care that the counter is "atomic", this is just a way to have a "reference-to-a-mutable-Int" to pass around
   */
  private def listAllGraphNodes(graph: Graph,
                                clusterCounter: AtomicInteger,
                                availableScatterVariables: Map[ScatterVariableNode, DotScatterVariableNode]): NodesAndLinks = {

    graph.nodes.toList foldMap {
      case ccn: CommandCallNode => NodesAndLinks(Set(DotCallNode(ccn)), upstreamLinks(ccn, DotCallNode(ccn), availableScatterVariables))
      case s: ScatterNode => handleScatter(s, clusterCounter, availableScatterVariables)
      case e: ExposedExpressionNode if hasCallAncestor(e) => NodesAndLinks(Set(DotExpressionNode(e)), upstreamLinks(e, DotExpressionNode(e), availableScatterVariables))
      case _ => nodeAndLinkMonoid.empty
    }


  }

  def handleScatter(scatterNode: ScatterNode,
                    clusterCounter: AtomicInteger,
                    knownScatterVariables: Map[ScatterVariableNode, DotScatterVariableNode]): NodesAndLinks = {
    val clusterNumber = clusterCounter.getAndIncrement()
    val id = s"cluster_$clusterNumber"

    def handleScatterVariableNode(scatterVariableNode: ScatterVariableNode): NodesAndLinks = {
      scatterNode.scatterVariableNodes foldMap { node =>
        val dotNode = DotScatterVariableNode.apply(node, clusterNumber)
        val links = upstreamLinks(node.linkToOuterGraph.graphNode, dotNode, Map.empty)
        NodesAndLinks(Set(dotNode), links)
      }
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
}

object GraphPrint {
  final case class DotLink(from: DotNode, to: DotNode) {
    def dotString = s"${from.id} -> ${to.id}"
  }

  final case class NodesAndLinks(nodes: Set[DotNode], links: Set[DotLink])

  final case class WorkflowDigraph(workflowName: String, digraph: NodesAndLinks) {
    def dotString: String =
      s"""|digraph $workflowName {
          |  rankdir=LR;
          |  compound=true;
          |
          |  # Links
          |  ${digraph.links.toList.flatMap(_.dotString.lines).mkString(System.lineSeparator + "  ")}
          |
          |  # Nodes
          |  ${digraph.nodes.toList.flatMap(_.dotString.lines).mkString(System.lineSeparator + "  ")}
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

  final case class DotExpressionNode(womType: WomType, valueName: String) extends DotNode {
    override def id: String = s"VALUE_$valueName"
    def dotString = s"""$id [shape="hexagon" label="${womType.toDisplayString} $valueName"]"""
  }
  object DotExpressionNode {
    def apply(expr: ExpressionNode): DotExpressionNode = DotExpressionNode(expr.womType, expr.identifier.localName.value)
  }

  final case class DotScatterVariableNode(womType: WomType, valueName: String, clusterNumber: Int) extends DotNode {
    override def id: String = s"SCATTER_${clusterNumber}_VARIABLE_$valueName"
    def dotString = s"""$id [shape="hexagon" label="scatter over ${womType.toDisplayString} as $valueName"]"""
  }
  object DotScatterVariableNode {
    def apply(svn: ScatterVariableNode, clusterNumber: Int): DotScatterVariableNode = DotScatterVariableNode(svn.womType, svn.identifier.localName.value, clusterNumber)
  }

  final case class DotScatterNode(id: String, nodes: Set[DotNode]) extends DotNode {
    override def dotString: String =
      s"""subgraph $id {
          |  rankdir=TB;
          |  style="filled,solid";
          |  fillcolor=white;
          |  ${nodes.toList.flatMap(_.dotString.lines).mkString(System.lineSeparator() + "  ")}
          |}""".stripMargin
  }

  case object CommandCallOutputPort {
    def unapply(port: OutputPort): Option[CommandCallNode] = port.graphNode match {
      case a: CommandCallNode => Some(a)
      case _ => None
    }
  }

  case object ExposedExpressionOutputPort {
    def unapply(port: OutputPort): Option[ExposedExpressionNode] = port.graphNode match {
      case a: ExposedExpressionNode => Some(a)
      case _ => None
    }
  }

  case object ExpressionOutputPort {
    def unapply(port: OutputPort): Option[ExpressionNode] = port.graphNode match {
      case a: ExpressionNode => Some(a)
      case _ => None
    }
  }

  case object OuterGraphInputPort {
    def unapply(port: OutputPort): Option[OuterGraphInputNode] = port.graphNode match {
      case a: OuterGraphInputNode => Some(a)
      case _ => None
    }
  }

  case object ScatterVariablePort {
    def unapply(port: OutputPort): Option[ScatterVariableNode] = port.graphNode match {
      case a: ScatterVariableNode => Some(a)
      case _ => None
    }
  }

  def upstreamLinks(originNode: GraphNode,
                    origin: DotNode,
                    availableScatterVariables: Map[ScatterVariableNode, DotScatterVariableNode]): Set[DotLink] = {
    def relevantAsUpstream(nodeToLink: GraphNode): Set[DotNode] = nodeToLink match {
      case ccn: CommandCallNode => Set(DotCallNode(ccn))
      case een: ExposedExpressionNode => if (hasCallAncestor(een)) { Set(DotExpressionNode(een)) } else Set.empty[DotNode]
      case en: ExpressionNode => upstreamLinksforNode(en)
      case svn: ScatterVariableNode => Set(availableScatterVariables(svn))
      case ogin: OuterGraphInputNode => relevantAsUpstream(ogin.linkToOuterGraph.graphNode)

      case _ => Set.empty[DotNode]
    }


    def upstreamLinksforNode(n: GraphNode) = n.upstreamPorts flatMap { upstreamPort =>
      println(s"Looking to make $upstreamPort a relevant upstream Node")
      upstreamPort match {
        case gatherPort: ScatterGathererPort =>
          relevantAsUpstream(gatherPort.outputToGather.singleUpstreamNode)
        case other => relevantAsUpstream(other.graphNode)
//        case CommandCallOutputPort(commandCallNode) =>
//          relevantAsUpstream(commandCallNode)
//        case ExposedExpressionOutputPort(exposedExpressionNode) =>
//          relevantAsUpstream(exposedExpressionNode)
//        case ExpressionOutputPort(expressionNode) =>
//          relevantAsUpstream(expressionNode)
//        case ScatterVariablePort(scatterVariableNode) =>
//          relevantAsUpstream(scatterVariableNode)
//        case OuterGraphInputPort(outerGraphInputNode) =>
//          relevantAsUpstream(outerGraphInputNode.linkToOuterGraph.graphNode)
//
//        case _ => Set.empty[DotNode]
      }
    }

    upstreamLinksforNode(originNode).map(DotLink(_, origin))
  }

  def hasCallAncestor(g: GraphNode) = g.upstreamAncestry.exists(_.isInstanceOf[CommandCallNode])
}
