package wom.views

import cats.implicits._
import cats.Monoid
import wom.executable.Executable
import wom.graph.expression.ExpressionNode
import wom.graph.{CommandCallNode, Graph, GraphNode}
import wom.views.GraphPrint._

final class GraphPrint(executableCallable: Executable) {

  val graph: Graph = executableCallable.graph

  def dotString = WorkflowDigraph(executableCallable.entryPoint.name, listAllGraphNodes).dotString

  implicit val nodeAndLinkMonoid: Monoid[NodesAndLinks] = cats.derived.MkMonoid[NodesAndLinks]

  private def listAllGraphNodes: NodesAndLinks = {


    def upstreamLinks(node: GraphNode, origin: CommandCallNode): Set[DotLink] = node.upstream flatMap {
      case commandCallNode: CommandCallNode => Set(DotLink(DotCallNode(commandCallNode), DotCallNode(origin)))
      case expression: ExpressionNode => upstreamLinks(expression, origin)
      case _ => Set.empty[DotLink]
    }

    graph.nodes.toList foldMap {
      case ccn: CommandCallNode => NodesAndLinks(Set(DotCallNode(ccn)), upstreamLinks(ccn, ccn))
      case _ => nodeAndLinkMonoid.empty
    }

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
          |  ${digraph.links.map(_.dotString).mkString(System.lineSeparator + "  ")}
          |  ${digraph.nodes.map(_.dotString).mkString(System.lineSeparator + "  ")}
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
}
