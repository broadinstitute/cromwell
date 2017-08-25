package wdltool.graph

import java.nio.file.{Files, Paths}
import java.util.UUID
import java.util.concurrent.atomic.AtomicInteger

import better.files.File
import cats.data.Validated.{Invalid, Valid}
import cats.derived.monoid._
import cats.derived.monoid.legacy._
import cats.implicits._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.cwl.CwlDecoder
import wdl4s.wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl4s.wom.executable.Executable
import wdl4s.wom.graph._
import wdltool.graph.WomGraph._

import scala.collection.JavaConverters._

class WomGraph(graphName: String, graph: Graph) {

  def indent(s: String) = s.lines.map(x => s"  $x").mkString(System.lineSeparator)
  def combine(ss: Iterable[String]) = ss.mkString(start="", sep=System.lineSeparator, end=System.lineSeparator)
  def indentAndCombine(ss: Iterable[String]) = combine(ss.map(indent))

  val digraphDot: String = {

    val workflowDigraph = WorkflowDigraph (dotSafe(graphName), listAllGraphNodes(graph))

    s"""|digraph ${workflowDigraph.workflowName}
        |{
        |  compound=true;
        |${indentAndCombine(workflowDigraph.digraph.links)}
        |${indentAndCombine(workflowDigraph.digraph.nodes)}
        |}
        |""".stripMargin

  }

  private lazy val clusterCount: AtomicInteger = new AtomicInteger(0)

  private[graph] def listAllGraphNodes(graph: Graph): NodesAndLinks = {
    graph.nodes foldMap nodesAndLinks
  }

  private def upstreamLinks(graphNode: GraphNode): Set[String] = for {
    inputPort <- graphNode.inputPorts
    upstreamPort = inputPort.upstream
  } yield s"${upstreamPort.graphId} -> ${inputPort.graphId}"

  private def portLine(p: GraphNodePort) = s"${p.graphId} [shape=${p.graphShape} label=${p.graphName}];"

  private def portLines(graphNode: GraphNode): String = graphNode match {
    case s: ScatterNode =>
      // Don't include the scatter expression input port here since they're added later in `internalScatterNodesAndLinks`
      // Round up the gathered output ports so it's obvious they're being gathered.
      s"""
        |${combine((s.inputPorts -- s.scatterVariableMapping.scatterInstantiatedExpression.inputPorts) map portLine)}
        |subgraph $nextCluster {
        |  style=filled;
        |  fillcolor=${s.graphFillColor}
        |  "${UUID.randomUUID}" [shape=plaintext label="gather ports"]
        |${indentAndCombine(s.outputPorts map portLine)}
        |}
      """.stripMargin
    case _ => combine((graphNode.outputPorts ++ graphNode.inputPorts) map portLine)
  }

  /**
    * For most graph nodes, we draw a single box for the node, containing the name and input and output ports.
    * For each input/output port, we include the links to upstream ports.
    */
  private def nodesAndLinks(graphNode: GraphNode): NodesAndLinks = {
    val internals = internalNodesAndLinks(graphNode)
    val clusterName = nextCluster

    val singleNode =
      s"""
         |subgraph $clusterName {
         |  style=filled;
         |  fillcolor=${graphNode.graphFillColor};
         |  ${graphNode.graphId} [shape=plaintext label=${graphNode.graphName}]
         |${indent(portLines(graphNode))}
         |${indentAndCombine(internals.nodes)}
         |}
         |""".stripMargin

    NodesAndLinks(Set(singleNode), internals.links ++ upstreamLinks(graphNode))
  }

  private def internalNodesAndLinks(graphNode: GraphNode): NodesAndLinks = graphNode match {
    case scatter: ScatterNode => internalScatterNodesAndLinks(scatter)
    case _ => NodesAndLinks.empty
  }

  def internalScatterNodesAndLinks(scatter: ScatterNode) = {
    val innerGraph = listAllGraphNodes(scatter.innerGraph) wrapNodes { n =>
      s"""
         |subgraph $nextCluster {
         |  style=filled;
         |  fillcolor=lightgray;
         |${indentAndCombine(n)}
         |}
         |""".stripMargin
    }
    val scatterVariableNodesAndLinks = {
      val scatterVariableExpressionId = s""""${UUID.randomUUID}""""
      val scatterVariableClusterName = nextCluster
      val node = s"""
         |subgraph $scatterVariableClusterName {
         |  style=filled;
         |  fillcolor=${scatter.graphFillColor};
         |  $scatterVariableExpressionId [shape=plaintext label="${scatter.scatterVariableMapping.graphInputNode.name} in ..."]
         |${indentAndCombine(scatter.scatterVariableMapping.scatterInstantiatedExpression.inputPorts map portLine)}
         |}
         """.stripMargin

      val scatterVariableInputLink = s"$scatterVariableExpressionId -> ${scatter.scatterVariableMapping.graphInputNode.singleOutputPort.graphId} [style=dashed ltail=$scatterVariableClusterName arrowhead=none]"
      val outputLinks = scatter.outputMapping map { outputPort => s"${outputPort.outputToGather.singleInputPort.graphId} -> ${outputPort.graphId} [style=dashed arrowhead=none]" }
      val links: Set[String] = outputLinks + scatterVariableInputLink

      NodesAndLinks(Set(node), links)
    }
    innerGraph |+| scatterVariableNodesAndLinks
  }

  private def nextCluster: String = "cluster_" + clusterCount.getAndIncrement()
}

object WomGraph {

  final case class WorkflowDigraph(workflowName: String, digraph: NodesAndLinks)
  final case class NodesAndLinks(nodes: Set[String], links: Set[String]) {
    def wrapNodes(f: Set[String] => String) = this.copy(Set(f(nodes)), links)
  }

  object NodesAndLinks {
    val empty = NodesAndLinks(Set.empty, Set.empty)
  }

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
    (for {
      clt <- CwlDecoder.decodeAllCwl(File(filePath)).
        value.
        unsafeRunSync
      wom <- clt.womExecutable.toEither
    } yield wom) match {
      case Right(womExecutable) => womExecutable
      case Left(e) => throw new Exception(s"Can't build WOM executable from CWL: ${e.toList.mkString("\n", "\n", "\n")}")
    }
  }
}
