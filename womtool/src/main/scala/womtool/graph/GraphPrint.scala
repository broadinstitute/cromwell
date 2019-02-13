package womtool.graph

import java.nio.file.{Files, Paths}
import java.util.concurrent.atomic.AtomicInteger

import cats.implicits._
import wdl.draft2.model.{Declaration, If, Scatter, WdlCall, WdlGraphNode, _}

import scala.collection.JavaConverters._

object GraphPrint {

  final case class WorkflowDigraph(workflowName: String, digraph: NodesAndLinks)
  final case class NodesAndLinks(nodes: Set[String], links: Set[String])
  implicit val monoid = cats.derived.MkMonoid[NodesAndLinks]

  def generateWorkflowDigraph(file: String): WorkflowDigraph = {
    // It's ok to use .get here, we're happy to throw an exception and crash the program!
    val namespace = WdlNamespaceWithWorkflow.load(Files.readAllLines(Paths.get(file)).asScala.mkString(System.lineSeparator()), Seq(WdlNamespace.fileResolver _)).get

    val digraph = listAllGraphNodes(namespace.workflow)

    WorkflowDigraph(namespace.workflow.unqualifiedName, digraph)
  }

  private val clusterCount: AtomicInteger = new AtomicInteger(0)


  private def listAllGraphNodes(scope: Scope): NodesAndLinks = {

    val callsAndDeclarations: Set[WdlGraphNode] = (scope.children collect {
      case w: WdlGraphNode if isCallOrCallBasedDeclaration(w) => w
    }).toSet

    val subGraphs: Set[WdlGraphNode] = (scope.children collect {
      case s: Scatter => s
      case i: If => i
    }).toSet

    def upstreamLinks(wdlGraphNode: WdlGraphNode, graphNodeName: String, suffix: String = ""): Set[String] = wdlGraphNode.upstream collect {
      case upstream: WdlGraphNode if isCallOrCallBasedDeclaration(upstream) =>
        val upstreamName = graphName(upstream)
        s""""$upstreamName" -> "$graphNodeName" $suffix"""
    }

    val thisLevelNodesAndLinks: NodesAndLinks = callsAndDeclarations.toList foldMap { graphNode =>
      val name = graphName(graphNode)
      val initialSet: Set[String] = graphNode match {
        case w: WdlGraphNode if isCallOrCallBasedDeclaration(w) => Set(s""""$name"""")
        case _ => Set.empty
      }

      NodesAndLinks(initialSet, upstreamLinks(graphNode, name))
    }

    val subGraphNodesAndLinks: NodesAndLinks = subGraphs.toList foldMap { wdlGraphNode =>
      val clusterName = "cluster_" + clusterCount.getAndIncrement()
      val subGraphName = graphName(wdlGraphNode)
      val subNodes = listAllGraphNodes(wdlGraphNode)
      val scope = s"""
         |subgraph $clusterName {
         |  ${subNodes.nodes.mkString(sep="\n  ")}
         |  "$subGraphName" [shape=plaintext]
         |}
      """.stripMargin

      NodesAndLinks(Set(scope), subNodes.links ++ upstreamLinks(wdlGraphNode, subGraphName, s"[lhead=$clusterName]"))
    }

    thisLevelNodesAndLinks |+| subGraphNodesAndLinks
  }

  private def isCallOrCallBasedDeclaration(w: WdlGraphNode): Boolean = w match {
    case _: WdlCall => true
    case w: Declaration if w.upstream.exists(isCallOrCallBasedDeclaration) => true
    case _ => false
  }

  private def dotSafe(s: String) = s.replaceAllLiterally("\"", "\\\"")

  private def graphName(g: WdlGraphNode): String = dotSafe(g match {
    case d: Declaration =>
      s"${d.womType.stableName} ${d.unqualifiedName}"
    case c: WdlCall =>
      s"call ${c.unqualifiedName}"
    case i: If =>
      s"if (${i.condition.toWomString})"
    case s: Scatter =>
      s"scatter (${s.collection.toWomString})"
    case c: CallOutput =>
      val exprString = c.expression.map(e => " = " + e.toWomString).getOrElse("")
      s"output { ${c.fullyQualifiedName}$exprString }"
    case other => s"${other.getClass.getSimpleName}: ${other.fullyQualifiedName}"
  })
}
