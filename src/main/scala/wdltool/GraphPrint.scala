package wdltool

import java.nio.file.{Files, Paths}
import java.util.concurrent.atomic.AtomicInteger

import wdl4s.wdl.{CallOutput, Declaration, If, Scatter, WdlCall, _}
import wdl4s.wdl.WdlGraphNode

import scala.collection.JavaConverters._

import cats.implicits._
import cats.derived.monoid._, legacy._

object GraphPrint {

  final case class WorkflowDigraph(workflowName: String, digraph: NodesAndLinks)
  final case class NodesAndLinks(nodes: Set[String], links: Set[String])

  def generateWorkflowDigraph(file: String, allNodesMode: Boolean): WorkflowDigraph = {
    // It's ok to use .get here, we're happy to throw an exception and crash the program!
    val namespace = WdlNamespaceWithWorkflow.load(Files.readAllLines(Paths.get(file)).asScala.mkString(System.lineSeparator()), Seq(WdlNamespace.fileResolver _)).get

    val digraph = listAllGraphNodes(namespace.workflow)

    WorkflowDigraph(namespace.workflow.unqualifiedName, digraph)
  }

  private val clusterCount: AtomicInteger = new AtomicInteger(0)


  private def listAllGraphNodes(scope: Scope): NodesAndLinks = {

    val callsAndDeclarations: Set[WdlGraphNode] = (scope.children collect {
      case w: WdlGraphNode if isCallOrDeclaration(w) => w
    }).toSet

    val subGraphs: Set[WdlGraphNode] = (scope.children collect {
      case s: Scatter => s
      case i: If => i
    }).toSet

    def upstreamLinks(wdlGraphNode: WdlGraphNode, graphNodeName: String, suffix: String = ""): Set[String] = wdlGraphNode.upstream collect {
      case upstream: WdlGraphNode if isCallOrDeclaration(upstream) =>
        val upstreamName = graphName(upstream)
        s""""$upstreamName" -> "$graphNodeName" $suffix"""
    }

    val thisLevelNodesAndLinks: NodesAndLinks = callsAndDeclarations foldMap { graphNode =>
      val name = graphName(graphNode)
      val initialSet: Set[String] = graphNode match {
        case w: WdlGraphNode if isCallOrDeclaration(w) => Set(s""""$name"""")
        case _ => Set.empty
      }

      val fromStart = if (graphNode.upstream.isEmpty) Set(s""""start" -> "$name"""") else Set.empty

      NodesAndLinks(initialSet, upstreamLinks(graphNode, name) ++ fromStart)
    }

    val subGraphNodesAndLinks: NodesAndLinks = subGraphs foldMap { wdlGraphNode =>
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

  private def isCallOrDeclaration(w: WdlGraphNode): Boolean = w.isInstanceOf[WdlCall] || w.isInstanceOf[Declaration]

  private def dotSafe(s: String) = s.replaceAllLiterally("\"", "\\\"")

  private def graphName(g: WdlGraphNode): String = dotSafe(g match {
    case d: Declaration =>
      val exprString = d.expression.map(e => " = " + e.toWdlString).getOrElse("")
      s"${d.wdlType.toWdlString} ${d.unqualifiedName}$exprString"
    case c: WdlCall =>
      s"call ${c.unqualifiedName}"
    case i: If =>
      s"if (${i.condition.toWdlString})"
    case s: Scatter =>
      s"scatter (${s.item} in ${s.collection.toWdlString})"
    case c: CallOutput =>
      val exprString = c.expression.map(e => " = " + e.toWdlString).getOrElse("")
      s"output { ${c.fullyQualifiedName}$exprString }"
    case other => s"${other.getClass.getSimpleName}: ${other.fullyQualifiedName}"
  })
}
