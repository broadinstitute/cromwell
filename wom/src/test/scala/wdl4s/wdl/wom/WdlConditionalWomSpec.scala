package wdl4s.wdl.wom

import cats.data.Validated.{Invalid, Valid}
import lenthall.collections.EnhancedCollections._
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.types._
import wdl4s.wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl4s.wom.graph._
import wdl4s.wom.graph.GraphNodePort.ConditionalOutputPort

class WdlConditionalWomSpec extends FlatSpec with Matchers {

  behavior of "WdlNamespaces with if blocks"

  it should "convert inputs from conditional tasks into workflow inputs, and conditional task outputs into workflow outputs" in {
    val conditionalTest =
      """task foo {
        |  Int i
        |  command { ... }
        |  output {
        |    String out = i
        |  }
        |}
        |
        |workflow conditional_test {
        |  Boolean b
        |  if (b) {
        |    call foo
        |  }
        |}""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(conditionalTest, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val conditionalTestGraph = namespace.womExecutable.flatMap(_.graph)

    conditionalTestGraph match {
      case Valid(g) => validateGraph(g)
      case Invalid(errors) => fail(s"Unable to build wom version of conditional foo from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateGraph(workflowGraph: Graph) = {

      case class OuterGraphValidations(conditionalNode: ConditionalNode, foo_i_inputNode: GraphInputNode)
      def validateOuterGraph: OuterGraphValidations = {
        val conditionalNode = workflowGraph.nodes.firstByType[ConditionalNode].getOrElse(fail("Resulting graph did not contain a ConditionalNode"))

        val inputNodes: Set[GraphInputNode] = workflowGraph.nodes.filterByType[GraphInputNode]

        val b_inputNode = inputNodes.find(_.name == "b").getOrElse(fail("Resulting graph did not contain the 'b' GraphInputNode"))
        b_inputNode.womType should be(WdlBooleanType)
        val foo_i_inputNode = inputNodes.find(_.name == "conditional_test.foo.i").getOrElse(fail("Resulting graph did not contain the 'foo.i' GraphInputNode"))
        foo_i_inputNode.womType should be(WdlIntegerType)

        val foo_out_output = workflowGraph.nodes.collectFirst {
          case gon: GraphOutputNode if gon.name == "foo.out" => gon
        }.getOrElse(fail("Resulting graph did not contain the 'foo.out' GraphOutputNode"))
        foo_out_output.womType should be(WdlOptionalType(WdlStringType))

        workflowGraph.nodes should be(Set(conditionalNode, foo_i_inputNode, b_inputNode, foo_out_output))
        OuterGraphValidations(conditionalNode, foo_i_inputNode)
      }

      case class InnerGraphValidations(foo_out_innerOutput: GraphOutputNode)
      def validateInnerGraph(validatedOuterGraph: OuterGraphValidations): InnerGraphValidations = {
        val foo_i_innerInput = validatedOuterGraph.conditionalNode.innerGraph.nodes.collectFirst {
          case gin: GraphInputNode if gin.name == "conditional_test.foo.i" => gin
        }.getOrElse(fail("Conditional inner graph did not contain a GraphInputNode 'foo.i'"))

        val foo_callNode = validatedOuterGraph.conditionalNode.innerGraph.nodes.collectFirst {
          case c: TaskCallNode if c.name == "foo" => c
        }.getOrElse(fail("Conditional inner graph did not contain a call to 'foo'"))

        val foo_out_innerOutput = validatedOuterGraph.conditionalNode.innerGraph.nodes.collectFirst {
          case gon: GraphOutputNode if gon.name == "foo.out" => gon
        }.getOrElse(fail("Conditional inner graph did not contain a GraphOutputNode 'foo.out'"))

        validatedOuterGraph.conditionalNode.innerGraph.nodes should be(Set(foo_i_innerInput, foo_callNode, foo_out_innerOutput))
        InnerGraphValidations(foo_out_innerOutput)
      }

      def validateConnections(validatedOuterGraph: OuterGraphValidations, validatedInnerGraph: InnerGraphValidations) = {
        // The ConditionalNode's output port is correctly associated with the inner graph's GraphOutputNode:
        validatedOuterGraph.conditionalNode.outputMapping.toList match {
          case ConditionalOutputPort(name, womType, outputToGather, _) :: Nil =>
            name should be("foo.out")
            womType should be(WdlOptionalType(WdlStringType))
            outputToGather eq validatedInnerGraph.foo_out_innerOutput should be(true)
          case other => fail("Expected exactly one output to be gathered in this conditional but got:" + other.mkString("\n", "\n", "\n"))
        }
      }

      val outer = validateOuterGraph
      val inner = validateInnerGraph(outer)
      validateConnections(outer, inner)
    }
  }


}
