package wdltool

import java.nio.file.{Files, Paths}

import wdl4s.{CallOutput, Declaration, If, Scatter, _}
import scala.collection.JavaConverters._

object GraphPrint {

  case class WorkflowDigraph(workflowName: String, digraph: Set[String])

  def generateWorkflowDigraph(file: String, allNodesMode: Boolean): WorkflowDigraph = {
    val namespace = WdlNamespaceWithWorkflow.load(Files.readAllLines(Paths.get(file)).asScala.mkString(System.lineSeparator()), Seq(WdlNamespace.fileResolver _))

    val digraph = if (allNodesMode) {
      listAllGraphNodes(namespace)
    } else {
      val executables = GraphPrint.listExecutableGraphNodes(namespace.workflow)
      listAllGraphNodes(namespace, graphNode => executables.contains(graphNode))
    }

    WorkflowDigraph(namespace.workflow.unqualifiedName, digraph)
  }

  private def defaultFilter: GraphNode => Boolean = _ => true

  private def listAllGraphNodes(namespace: WdlNamespaceWithWorkflow, filter: GraphNode => Boolean = defaultFilter): Set[String] = {

    val graphNodes = namespace.descendants collect {
      case g: GraphNode if filter(g) => g
    }

    graphNodes flatMap { graphNode =>
      val name = graphName(graphNode)
      val initialSet: Set[String] = graphNode match {
        case c: Call => Set(s""""${dotSafe(name)}"""")
        case _ => Set.empty
      }
      val upstreamLinks = graphNode.upstream collect {
        case upstream if filter(upstream) =>
          val upstreamName = graphName(upstream)
          s""""${dotSafe(upstreamName)}" -> "${dotSafe(name)}""""
      }

      initialSet ++ upstreamLinks
    }
  }

  private def listExecutableGraphNodes(s: Scope): Set[GraphNode] = {
    s.children.toSet flatMap { child: Scope => child match {
      case call: Call => Set[GraphNode](call)
      case scatter: Scatter => Set[GraphNode](scatter) ++ listExecutableGraphNodes(scatter)
      case i: If => Set[GraphNode](i) ++ listExecutableGraphNodes(i)
      case declaration: Declaration => Set[GraphNode](declaration)
      case _ => Set.empty[GraphNode]
    }}
  }


  private def dotSafe(s: String) = s.replaceAllLiterally("\"", "\\\"")

  private def graphName(g: GraphNode): String = g match {
    case d: Declaration =>
      val exprString = d.expression.map(e => " = " + e.toWdlString).getOrElse("")
      s"${d.wdlType.toWdlString} ${d.fullyQualifiedName}$exprString"
    case c: Call =>
      s"call ${c.fullyQualifiedName}"
    case i: If =>
      s"if (${i.condition.toWdlString})"
    case s: Scatter =>
      s"scatter (${s.item} in ${s.collection.toWdlString})"
    case c: CallOutput =>
      val exprString = c.expression.map(e => " = " + e.toWdlString).getOrElse("")
      s"output { ${c.fullyQualifiedName}$exprString }"
    case other => s"${other.getClass.getSimpleName}: ${other.fullyQualifiedName}"
  }
}
