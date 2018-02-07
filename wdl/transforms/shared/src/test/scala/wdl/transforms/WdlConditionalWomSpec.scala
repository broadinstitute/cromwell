package wdl.transforms

import cats.data.Validated.{Invalid, Valid}
import common.collections.EnhancedCollections._
import org.scalatest.{FlatSpec, Matchers}
import wdl.{WdlNamespace, WdlNamespaceWithWorkflow, WdlWorkflow}
import wom.transforms.WomWorkflowDefinitionMaker
import wom.graph.GraphNodePort.ConditionalOutputPort
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.types.{WomBooleanType, WomIntegerType, WomOptionalType, WomStringType}
import wom.transforms.WomWorkflowDefinitionMaker.ops._

class WdlConditionalWomSpec(implicit workflowDefinitionMaker: WomWorkflowDefinitionMaker[WdlWorkflow]) extends FlatSpec with Matchers {

  behavior of "WdlNamespaces with if blocks"

  it should "respect expression inputs and convert conditional task outputs into workflow outputs" in {
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
    val conditionalTestGraph = namespace.workflow.toWomWorkflowDefinition.map(_.graph)

    conditionalTestGraph match {
      case Valid(g) => validateGraph(g)
      case Invalid(errors) => fail(s"Unable to build wom version of conditional foo from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateGraph(workflowGraph: Graph) = {

      case class OuterGraphValidations(conditionalNode: ConditionalNode, foo_i_inputNode: GraphInputNode)
      def validateOuterGraph: OuterGraphValidations = {
        val conditionalNode = workflowGraph.nodes.firstByType[ConditionalNode].getOrElse(fail("Resulting graph did not contain a ConditionalNode"))

        val inputNodes: Set[GraphInputNode] = workflowGraph.nodes.filterByType[GraphInputNode]

        val b_inputNode = inputNodes.find(_.localName == "b").getOrElse(fail("Resulting graph did not contain the 'b' GraphInputNode"))
        b_inputNode.womType should be(WomBooleanType)
        val foo_i_inputNode = inputNodes.find(_.localName == "foo.i").getOrElse(fail("Resulting graph did not contain the 'foo.i' GraphInputNode"))
        foo_i_inputNode.womType should be(WomIntegerType)

        val foo_out_output = workflowGraph.nodes.collectFirst {
          case gon: GraphOutputNode if gon.localName == "foo.out" => gon
        }.getOrElse(fail("Resulting graph did not contain the 'foo.out' GraphOutputNode"))
        foo_out_output.womType should be(WomOptionalType(WomStringType))
        foo_out_output.identifier.fullyQualifiedName.value shouldBe "conditional_test.foo.out"
        
        val expressionNode = workflowGraph.nodes.collectFirst {
          case expr: ExpressionNode if expr.localName == "conditional" => expr
        }.getOrElse(fail("Resulting graph did not contain the 'conditional' ExpressionNode"))

        workflowGraph.nodes should be(Set(conditionalNode, foo_i_inputNode, b_inputNode, foo_out_output, expressionNode))
        OuterGraphValidations(conditionalNode, foo_i_inputNode)
      }

      case class InnerGraphValidations(foo_out_innerOutput: GraphOutputNode)
      def validateInnerGraph(validatedOuterGraph: OuterGraphValidations): InnerGraphValidations = {
        val foo_i_innerInput = validatedOuterGraph.conditionalNode.innerGraph.nodes.collectFirst {
          case gin: ExternalGraphInputNode if gin.identifier.fullyQualifiedName.value == "conditional_test.foo.i" => gin
        }.getOrElse(fail("Conditional inner graph did not contain a GraphInputNode 'foo.i'"))

        val foo_callNode = validatedOuterGraph.conditionalNode.innerGraph.nodes.collectFirst {
          case c: CommandCallNode if c.localName == "foo" => c
        }.getOrElse(fail("Conditional inner graph did not contain a call to 'foo'"))

        foo_callNode.identifier.fullyQualifiedName.value shouldBe "conditional_test.foo"

        val foo_out_innerOutput = validatedOuterGraph.conditionalNode.innerGraph.nodes.collectFirst {
          case gon: GraphOutputNode if gon.localName == "foo.out" => gon
        }.getOrElse(fail("Conditional inner graph did not contain a GraphOutputNode 'foo.out'"))

        validatedOuterGraph.conditionalNode.innerGraph.nodes should be(Set(foo_i_innerInput, foo_callNode, foo_out_innerOutput))
        InnerGraphValidations(foo_out_innerOutput)
      }

      def validateConnections(validatedOuterGraph: OuterGraphValidations, validatedInnerGraph: InnerGraphValidations) = {
        // The ConditionalNode's output port is correctly associated with the inner graph's GraphOutputNode:
        validatedOuterGraph.conditionalNode.conditionalOutputPorts.toList match {
          case (port @ ConditionalOutputPort(outputToGather, _)) :: Nil =>
            port.name should be("foo.out")
            port.womType should be(WomOptionalType(WomStringType))
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
