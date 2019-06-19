package wom.graph

import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ShortCircuitingFlatMap
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wom.RuntimeAttributes
import wom.callable.Callable.{OutputDefinition, RequiredInputDefinition}
import wom.callable.CallableTaskDefinition
import wom.expression.PlaceholderWomExpression
import wom.graph.CallNode.{CallNodeAndNewNodes, CallNodeBuilder, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.{AnonymousExpressionNode, PlainAnonymousExpressionNode}
import wom.types.{WomArrayType, WomIntegerType, WomStringType}

class ScatterNodeSpec extends FlatSpec with Matchers {
  behavior of "ScatterNode"

  val fooInputDef = RequiredInputDefinition("i", WomIntegerType)
  val task_foo = CallableTaskDefinition(
    name = "foo",
    commandTemplateBuilder = null,
    runtimeAttributes = RuntimeAttributes(Map.empty),
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List(OutputDefinition("out", WomStringType, PlaceholderWomExpression(Set.empty, WomStringType))),
    inputs = List(fooInputDef),
    adHocFileCreation = Set.empty,
    environmentExpressions = Map.empty,
    sourceLocation = None
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
    val xs_inputNode = RequiredGraphInputNode(WomIdentifier("xs"), WomArrayType(WomIntegerType), "xs")

    val xsExpression = PlaceholderWomExpression(Set("xs"), WomArrayType(WomIntegerType))
    val xsExpressionAsInput = AnonymousExpressionNode
      .fromInputMapping(WomIdentifier("x"), xsExpression, Map("xs" -> xs_inputNode.singleOutputPort), PlainAnonymousExpressionNode.apply)
      .valueOr(failures => fail(s"Failed to create expression node: ${failures.toList.mkString(", ")}"))

    val x_inputNode = ScatterVariableNode(WomIdentifier("x"), xsExpressionAsInput, WomArrayType(WomIntegerType))
    val fooNodeBuilder = new CallNodeBuilder()
    val fooInputFold = InputDefinitionFold(
      mappings = List(
        fooInputDef -> Coproduct[InputDefinitionPointer](x_inputNode.singleOutputPort: OutputPort)
      ),
      callInputPorts = Set(
        fooNodeBuilder.makeInputPort(fooInputDef, x_inputNode.singleOutputPort)
      ),
      newGraphInputNodes = Set.empty
    )
    val CallNodeAndNewNodes(foo_callNode, _, _, _) = fooNodeBuilder.build(WomIdentifier("foo"), task_foo, fooInputFold, Set.empty, None)
    val foo_call_outNode = PortBasedGraphOutputNode(WomIdentifier("foo.out"), WomStringType, foo_callNode.outputByName("foo.out").getOrElse(fail("foo CallNode didn't contain the expected 'out' output")))
    val scatterGraph = Graph.validateAndConstruct(Set(foo_callNode, x_inputNode, foo_call_outNode)) match {
      case Valid(sg) => sg
      case Invalid(es) => fail("Failed to make scatter graph: " + es.toList.mkString(", "))
    }

    val scatterNodeWithInputs = ScatterNode.scatterOverGraph(
      scatterGraph,
      x_inputNode
    )

    val scatterNode = scatterNodeWithInputs.node
    scatterNodeWithInputs.newInputs.size should be(0)

    val workflowGraphValidation = for {
      foo_scatterOutput <- scatterNode.outputByName("foo.out")
      z_workflowOutput = PortBasedGraphOutputNode(WomIdentifier("z"), WomArrayType(WomStringType), foo_scatterOutput)
      graph <- Graph.validateAndConstruct(scatterNode.nodes ++ Set(xs_inputNode, z_workflowOutput))
    } yield graph

    workflowGraphValidation match {
      case Valid(wg) => validate(wg, scatterNode)
      case Invalid(es) => fail("Failed to make workflow graph: " + es.toList.mkString(", "))
    }

    /**
      * WDL -> WOM conversion has succeeded, now let's check we built it right!
      */
    def validate(workflowGraph: Graph, scatterNode: ScatterNode) = {
      workflowGraph.nodes collect { case gin: GraphInputNode => gin.localName } should be(Set("xs"))
      workflowGraph.nodes collect { case gon: PortBasedGraphOutputNode => gon.localName } should be(Set("z"))
      workflowGraph.nodes collect { case cn: CallNode => cn.localName } should be(Set.empty)
      (workflowGraph.nodes collect { case sn: ScatterNode => sn }).size should be(1)

      // The output links back to the scatter:
      val finalOutput = (workflowGraph.nodes collect { case gon: PortBasedGraphOutputNode => gon }).head
      finalOutput.upstream should be(Set(scatterNode))

      // The scatter output port is called foo.out:
      finalOutput.inputPorts.head.upstream.name should be("foo.out")

      // foo.out links back to the correct output in the inner graph:
      val innerGraphFooOutNode = scatterNode.outputMapping.find(_.name == "foo.out").getOrElse(fail("Scatter couldn't link back the foo.out output.")).outputToGather
      innerGraphFooOutNode.womType should be(WomStringType)
      innerGraphFooOutNode.upstream.size should be(1)
      innerGraphFooOutNode.upstream.head match {
        case c: CallNode => c.localName should be("foo")
        case _ => fail("Expected the inner graph to have a Call Node!")
      }
    }
  }
}
