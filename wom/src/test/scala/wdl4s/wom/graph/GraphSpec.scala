package wdl4s.wom.graph

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.RuntimeAttributes
import wdl4s.wdl.types.{WdlFileType, WdlIntegerType, WdlStringType}
import wdl4s.wom.callable.Callable.{OutputDefinition, RequiredInputDefinition}
import wdl4s.wom.callable.{TaskDefinition, WorkflowDefinition}
import wdl4s.wom.graph.CallNode.CallWithInputs

class GraphSpec extends FlatSpec with Matchers {
  behavior of "Graph"

  def makeThreeStep: Graph = {
    val taskDefinition_ps = TaskDefinition(
      name = "ps",
      commandTemplate = null,
      runtimeAttributes = new RuntimeAttributes(attrs = Map.empty),
      meta = Map.empty,
      parameterMeta = Map.empty,
      outputs = Set(OutputDefinition("procs", WdlFileType, null)),
      inputs = Set.empty,
      declarations = List.empty
    )

    val taskDefinition_cgrep = TaskDefinition(
      name = "cgrep",
      commandTemplate = null,
      runtimeAttributes = new RuntimeAttributes(attrs = Map.empty),
      meta = Map.empty,
      parameterMeta = Map.empty,
      outputs = Set(OutputDefinition("count", WdlIntegerType, null)),
      inputs = Set(RequiredInputDefinition("pattern", WdlStringType), RequiredInputDefinition("in_file", WdlFileType)),
      declarations = List.empty
    )

    val taskDefinition_wc = TaskDefinition(
      name = "wc",
      commandTemplate = null,
      runtimeAttributes = new RuntimeAttributes(attrs = Map.empty),
      meta = Map.empty,
      parameterMeta = Map.empty,
      outputs = Set(OutputDefinition("count", WdlIntegerType, null)),
      inputs = Set(RequiredInputDefinition("in_file", WdlFileType)),
      declarations = List.empty
    )

    val CallWithInputs(psCall, psGraphInputs) = CallNode.callWithInputs("ps", taskDefinition_ps, Map.empty)
    val ps_procsOutputPort = psCall.outputByName("procs").getOrElse(fail("Unexpectedly unable to find 'procs' output"))

    val CallWithInputs(cgrepCall, cgrepGraphInputs) = CallNode.callWithInputs("cgrep", taskDefinition_cgrep, Map("in_file" -> ps_procsOutputPort))
    val CallWithInputs(wcCall, wcGraphInputs) = CallNode.callWithInputs("wc", taskDefinition_wc, Map("in_file" -> ps_procsOutputPort))

    val graphNodes: Set[GraphNode] =
      Set[GraphNode](psCall, cgrepCall, wcCall)
        .union(psGraphInputs.toSet[GraphNode])
        .union(cgrepGraphInputs.toSet[GraphNode])
        .union(wcGraphInputs.toSet[GraphNode])

    Graph.validateAndConstruct(graphNodes) match {
      case Valid(wg) => wg.withDefaultOutputs
      case Invalid(errors) => fail(s"Unable to validate graph: ${errors.toList.mkString("\n", "\n", "\n")}")
    }
  }

  it should "be able to represent three step" in {
    val workflowGraph = makeThreeStep

    workflowGraph.nodes collect { case gin: GraphInputNode => gin.name } should be(Set("cgrep.pattern"))
    workflowGraph.nodes collect { case gon: GraphOutputNode => gon.name } should be(Set("wc.count", "cgrep.count", "ps.procs"))
    workflowGraph.nodes collect { case cn: CallNode => cn.name } should be(Set("wc", "cgrep", "ps"))
  }

  it should "be able to represent calls to sub-workflows" in {
    val threeStepWorkflow = WorkflowDefinition("three_step", makeThreeStep, Map.empty, Map.empty, List.empty)
    val CallWithInputs(threeStepCall, threeStepInputs) = CallNode.callWithInputs("three_step", threeStepWorkflow, Map.empty)

    val workflowGraph = Graph.validateAndConstruct(Set[GraphNode](threeStepCall).union(threeStepInputs.toSet[GraphNode])) match {
      case Valid(wg) => wg.withDefaultOutputs
      case Invalid(errors) => fail(s"Unable to validate graph: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    workflowGraph.nodes collect { case gin: GraphInputNode => gin.name } should be(Set("three_step.cgrep.pattern"))
    workflowGraph.nodes collect { case gon: GraphOutputNode => gon.name } should be(Set("three_step.wc.count", "three_step.cgrep.count", "three_step.ps.procs"))
    workflowGraph.nodes collect { case cn: CallNode => cn.name } should be(Set("three_step"))
  }
}
