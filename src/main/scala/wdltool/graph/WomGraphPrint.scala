package wdltool.graph

import java.nio.file.{Files, Paths}
import java.util.concurrent.atomic.AtomicInteger

import cats.data.Validated.{Invalid, Valid}
import cats.derived.monoid._
import cats.derived.monoid.legacy._
import cats.implicits._
import wdl4s.cwl.{CommandLineTool, CwlCodecs, Workflow, Yaml}
import wdl4s.wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl4s.wom.executable.Executable
import wdl4s.wom.graph.{Graph, GraphNode}

import scala.collection.JavaConverters._
import wdltool.graph._

object WomGraphPrint {

  final case class WorkflowDigraph(workflowName: String, digraph: NodesAndLinks)
  final case class NodesAndLinks(nodes: Set[String], links: Set[String])

  def readFile(filePath: String): String = Files.readAllLines(Paths.get(filePath)).asScala.mkString(System.lineSeparator())

  def womExecutableFromWdl(filePath: String): Executable = {
    val namespace = WdlNamespaceWithWorkflow.load(readFile(filePath), Seq(WdlNamespace.fileResolver _)).get
    namespace.womExecutable
  }

  def womExecutableFromCwl(filePath: String): Executable = {
    val yaml: Yaml = readFile(filePath)
    CwlCodecs.decodeCwl(yaml) match {
      case Right((clt: CommandLineTool, _)) =>
        clt.womExecutable match {
          case Valid(wom) => wom
          case Invalid(e) => throw new Exception(s"Can't build WOM executable from CWL CommandLineTool: ${e.toList.mkString("\n", "\n", "\n")}")
        }
      case Right((wf: Workflow, m)) if m.isEmpty =>
        wf.womExecutable(Map.empty) match {
          case Valid(wom) => wom
          case Invalid(e) => throw new Exception(s"Can't build WOM executable from CWL Workflow: ${e.toList.mkString("\n", "\n", "\n")}")
        }
        // TODO: Don't throw away the error information! Might as well return it:
      case other => throw new Exception(s"Can't read $filePath")
    }
  }

  def generateWorkflowDigraph(file: String, auxFiles: Seq[String]): WorkflowDigraph = {

    val womExecutable = if (file.toLowerCase().endsWith("wdl")) womExecutableFromWdl(file) else womExecutableFromCwl(file)
    womExecutable.graph match {
      case Valid(graph) =>
        val digraph = listAllGraphNodes (graph)
        WorkflowDigraph (dotSafe(womExecutable.entryPoint.name), digraph)
      case Invalid(errors) =>
        throw new IllegalArgumentException(errors.toList.mkString("\n", "\n", "\n"))
    }
  }

  private val clusterCount: AtomicInteger = new AtomicInteger(0)

  private def listAllGraphNodes(graph: Graph): NodesAndLinks = {

    graph.nodes foldMap nodesAndLinks
  }

  private def nodesAndLinks(graphNode: GraphNode): NodesAndLinks = graphNode match {
    case other => defaultNodesAndLinks(graphNode)
  }

  /**
    * For most graph nodes, we draw a single box for the node, containing the name and input and output ports.
    * For each input/output port, we include the links to upstream ports.
    */
  private def defaultNodesAndLinks(graphNode: GraphNode): NodesAndLinks = {
    val clusterName = "cluster_" + clusterCount.getAndIncrement()
    val portNodes = (graphNode.outputPorts ++ graphNode.inputPorts) map { p =>
      s"${p.graphId} [shape=${p.graphShape} label=${p.graphName}];"
    }
    val links = for {
      inputPort <- graphNode.inputPorts
      upstreamPort = inputPort.upstream
    } yield s"${upstreamPort.graphId} -> ${inputPort.graphId}"

    val singleNode =
      s"""
         |subgraph $clusterName {
         |  style=filled;
         |  fillcolor=${graphNode.graphFillColor};
         |  ${portNodes.mkString(sep="\n  ")}
         |  ${graphNode.graphId} [shape=plaintext label=${graphNode.graphName}]
         |}
         |""".stripMargin

    NodesAndLinks(Set(singleNode), links)
  }
}
