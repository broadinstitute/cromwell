package wdltool.graph

import java.nio.file.{Files, Paths}
import java.util.concurrent.atomic.AtomicInteger

import cats.data.Validated.{Invalid, Valid}
import cats.derived.monoid._
import cats.derived.monoid.legacy._
import cats.implicits._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.cwl.{CommandLineTool, CwlCodecs, Workflow}
import wdl4s.wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl4s.wom.executable.Executable
import wdl4s.wom.graph.{Graph, GraphNode}
import wdltool.graph.WomGraph._

import scala.collection.JavaConverters._

class WomGraph(graphName: String, graph: Graph) {

  val digraphDot: String = {

    val workflowDigraph = WorkflowDigraph (dotSafe(graphName), listAllGraphNodes(graph))

    s"""|digraph ${workflowDigraph.workflowName}
        |{
        |compound=true;
        |${workflowDigraph.digraph.links.mkString(System.lineSeparator + "  ")}
        |${workflowDigraph.digraph.nodes.mkString(System.lineSeparator + "  ")}
        |}
        |""".stripMargin

  }

  private lazy val clusterCount: AtomicInteger = new AtomicInteger(0)

  private[graph] def listAllGraphNodes(graph: Graph): NodesAndLinks = {

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

object WomGraph {

  final case class WorkflowDigraph(workflowName: String, digraph: NodesAndLinks)
  final case class NodesAndLinks(nodes: Set[String], links: Set[String])

  def fromFiles(mainFile: String, auxFiles: Seq[String]): ErrorOr[WomGraph] = {
    val womExecutable = if (mainFile.toLowerCase().endsWith("wdl")) womExecutableFromWdl(mainFile) else womExecutableFromCwl(mainFile)
    womExecutable.graph map { graph => new WomGraph(womExecutable.entryPoint.name, graph) }
  }

  private def readFile(filePath: String): String = Files.readAllLines(Paths.get(filePath)).asScala.mkString(System.lineSeparator())

  private def womExecutableFromWdl(filePath: String): Executable = {
    val namespace = WdlNamespaceWithWorkflow.load(readFile(filePath), Seq(WdlNamespace.fileResolver _)).get
    namespace.womExecutable match {
      case Valid(wom) => wom
      case Invalid(e) => throw new Exception(s"Can't build WOM executable from WDL namespace: ${e.toList.mkString("\n", "\n", "\n")}")
    }
  }

  private def womExecutableFromCwl(filePath: String): Executable = {
    val yaml: String = readFile(filePath)
    CwlCodecs.decodeCwl(yaml) match {
      case Valid((clt: CommandLineTool, _)) =>
        clt.womExecutable match {
          case Valid(wom) => wom
          case Invalid(e) => throw new Exception(s"Can't build WOM executable from CWL CommandLineTool: ${e.toList.mkString("\n", "\n", "\n")}")
        }
      case Valid((wf: Workflow, m)) if m.isEmpty =>
        wf.womExecutable(Map.empty) match {
          case Valid(wom) => wom
          case Invalid(e) => throw new Exception(s"Can't build WOM executable from CWL Workflow: ${e.toList.mkString("\n", "\n", "\n")}")
        }
      case Invalid(errors) => throw new Exception(s"Can't read $filePath")
    }
  }
}
