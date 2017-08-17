package wdl4s.wom.graph

import cats.data.Validated.{Invalid, Valid}
import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.RuntimeAttributes
import wdl4s.wdl.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.wom.callable.Callable.{OutputDefinition, RequiredInputDefinition}
import wdl4s.wom.callable.TaskDefinition
import wdl4s.wom.expression.PlaceholderWomExpression
import wdl4s.wom.graph.CallNode.CallWithInputs
import wdl4s.wom.graph.ScatterNode.ScatterNodeWithInputs

class ScatterNodeSpec extends FlatSpec with Matchers {
  behavior of "ScatterNode"

  val task_foo = TaskDefinition(name = "foo",
    commandTemplate = null,
    runtimeAttributes = new RuntimeAttributes(Map.empty),
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = Set(OutputDefinition("out", WdlStringType, PlaceholderWomExpression(Set.empty, WdlStringType))),
    inputs = Set(RequiredInputDefinition("i", WdlIntegerType)),
    declarations = List.empty
  )

  /**
    * Representing, approximately:
    *
    * task foo {
    *   Int i
    *   command { ... }
    *   output {
    *     String out = ...
    *   }
    *
    * workflow {
    *   Array[Int] xs
    *   scatter (x in xs) {
    *     call foo { input: i = x }
    *   }
    *   output {
    *     Array[String] z = foo.out
    *   }
    * }
    *
    */
  it should "be able to wrap a single task call" in {
    val xs_inputNode = RequiredGraphInputNode("xs", WdlArrayType(WdlIntegerType))

    val x_inputNode = RequiredGraphInputNode("x", WdlIntegerType)
    val CallWithInputs(foo_callNode, _) = CallNode.callWithInputs("foo", task_foo, Map("i" -> x_inputNode.singleOutputPort), Set.empty).getOrElse(fail("Unable to call foo_callNode"))
    val scatterGraph = Graph.validateAndConstruct(Set(foo_callNode, x_inputNode)) match {
      case Valid(sg) => sg.withDefaultOutputs
      case Invalid(es) => fail("Failed to make scatter graph: " + es.toList.mkString(", "))
    }

    val ScatterNodeWithInputs(scatterNode, unsuppliedScatterInputNodes) = ScatterNode.scatterOverGraph(scatterGraph, xs_inputNode.singleOutputPort, x_inputNode, Map.empty)
    unsuppliedScatterInputNodes.size should be(0)

    val workflowGraphValidation = for {
      foo_scatterOutput <- scatterNode.outputByName("foo.out")
      z_workflowOutput = PortBasedGraphOutputNode("z", WdlArrayType(WdlStringType), foo_scatterOutput)
      graph <- Graph.validateAndConstruct(Set(scatterNode, xs_inputNode, z_workflowOutput))
    } yield graph

    val workflowGraph = workflowGraphValidation match {
      case Valid(wg) => wg
      case Invalid(es) => fail("Failed to make workflow graph: " + es.toList.mkString(", "))
    }

    //
    // Above: Tests that we can create the scatter
    // Below: Tests that it was made right!
    //


    workflowGraph.nodes collect { case gin: GraphInputNode => gin.name } should be(Set("xs"))
    workflowGraph.nodes collect { case gon: PortBasedGraphOutputNode => gon.name } should be(Set("z"))
    workflowGraph.nodes collect { case cn: CallNode => cn.name } should be(Set.empty)
    (workflowGraph.nodes collect { case sn: ScatterNode => sn }).size should be(1)

    // The output links back to the scatter:
    val finalOutput = (workflowGraph.nodes collect { case gon: PortBasedGraphOutputNode => gon }).head
    finalOutput.upstream should be(Set(scatterNode))

    // The scatter output port is called foo.out:
    finalOutput.inputPorts.head.upstream.name should be("foo.out")

    // foo.out links back to the correct output in the inner graph:
    val innerGraphFooOutNode = scatterNode.outputMapping.find(_.name == "foo.out").getOrElse(fail("Scatter couldn't link back the foo.out output.")).outputToGather
    innerGraphFooOutNode.womType should be(WdlStringType)
    innerGraphFooOutNode.upstream.size should be(1)
    innerGraphFooOutNode.upstream.head match {
      case c: CallNode => c.name should be("foo")
      case _ => fail("Expected the inner graph to have a Call Node!")
    }

  }
}
