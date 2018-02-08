package wdl.transforms.wdlwom

import cats.data.Validated.{Invalid, Valid}
import common.collections.EnhancedCollections._
import org.scalatest.{FlatSpec, Matchers}
import wom.graph.GraphNodePort.ScatterGathererPort
import wom.graph.expression.ExpressionNode
import wom.graph.{GraphInputNode, ScatterNode, _}
import wom.transforms.WomWorkflowDefinitionMaker
import wom.types.{WomArrayType, WomIntegerType, WomStringType}

class WdlScatterWomSpec(implicit workflowDefinitionMaker: WomWorkflowDefinitionMaker[WdlWorkflow]) extends FlatSpec with Matchers {

  behavior of "WdlNamespaces with scatters"

  it should "convert gathered scatter outputs to workflow outputs" in {
    val scatterTest =
      """task foo {
        |  Int i
        |  command { ... }
        |  output {
        |    String out = i
        |  }
        |}
        |
        |workflow scatter_test {
        |  Array[Int] xs
        |  scatter (x in xs) {
        |    call foo { input: i = x }
        |  }
        |}""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(scatterTest, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    val scatterTestGraph = namespace.workflow.toWomWorkflowDefinition.map(_.graph)

    scatterTestGraph match {
      case Valid(g) => validateGraph(g)
      case Invalid(errors) => fail(s"Unable to build wom version of scatter foo from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateGraph(workflowGraph: Graph) = {

      case class OuterGraphValidations(scatterNode: ScatterNode, xs_inputNode: GraphInputNode)
      def validateOuterGraph: OuterGraphValidations = {
        val scatterNode = workflowGraph.nodes.firstByType[ScatterNode].getOrElse(fail("Resulting graph did not contain a ScatterNode"))

        val xs_inputNode = workflowGraph.nodes.collectFirst {
          case gin: GraphInputNode if gin.localName == "xs" => gin
        }.getOrElse(fail("Resulting graph did not contain the 'xs' GraphInputNode"))

        val scatterExpressionNode = workflowGraph.nodes.collectFirst {
          case expr: ExpressionNode if expr.localName == "x" => expr
        }.getOrElse(fail("Resulting graph did not contain the 'x' ExpressionNode"))
        
        scatterNode.inputPorts.map(_.upstream) shouldBe Set(scatterExpressionNode.singleExpressionOutputPort)
        
        val foo_out_output = workflowGraph.nodes.collectFirst {
          case gon: GraphOutputNode if gon.localName == "foo.out" => gon
        }.getOrElse(fail("Resulting graph did not contain the 'foo.out' GraphOutputNode"))
        foo_out_output.womType should be(WomArrayType(WomStringType))
        foo_out_output.identifier.fullyQualifiedName.value shouldBe "scatter_test.foo.out"

        workflowGraph.nodes should be(Set(scatterNode, xs_inputNode, foo_out_output, scatterExpressionNode))
        OuterGraphValidations(scatterNode, xs_inputNode)
      }

      case class InnerGraphValidations(x_scatterCollectionInput: GraphInputNode, foo_out_innerOutput: GraphOutputNode)
      def validateInnerGraph(validatedOuterGraph: OuterGraphValidations): InnerGraphValidations = {
        val x_scatterCollectionInput = validatedOuterGraph.scatterNode.innerGraph.nodes.collectFirst {
          case gin: GraphInputNode if gin.localName == "x" => gin
        }.getOrElse(fail("Scatter inner graph did not contain a GraphInputNode 'x'"))

        val foo_callNode = validatedOuterGraph.scatterNode.innerGraph.nodes.collectFirst {
          case c: CommandCallNode if c.localName == "foo" => c
        }.getOrElse(fail("Scatter inner graph did not contain a call to 'foo'"))
        
        foo_callNode.identifier.fullyQualifiedName.value shouldBe "scatter_test.foo"

        val foo_out_innerOutput = validatedOuterGraph.scatterNode.innerGraph.nodes.collectFirst {
          case gon: GraphOutputNode if gon.localName == "foo.out" => gon
        }.getOrElse(fail("Scatter inner graph did not contain a GraphOutputNode 'foo.out'"))

        val foo_out_i_expressionNode = validatedOuterGraph.scatterNode.innerGraph.nodes.collectFirst {
          case expr: ExpressionNode if expr.localName == "foo.i" => expr
        }.getOrElse(fail("Scatter inner graph did not contain a ExpressionNode 'scatter_test.foo.i'"))

        validatedOuterGraph.scatterNode.innerGraph.nodes should be(Set(x_scatterCollectionInput, foo_callNode, foo_out_innerOutput, foo_out_i_expressionNode))
        InnerGraphValidations(x_scatterCollectionInput, foo_out_innerOutput)
      }

      def validateConnections(validatedOuterGraph: OuterGraphValidations, validatedInnerGraph: InnerGraphValidations) = {
        // The scatter collection links to its predecessor
        validatedOuterGraph.scatterNode.scatterCollectionExpressionNodes.head.inputPorts.map(_.upstream.graphNode) should be(Set(validatedOuterGraph.xs_inputNode))

        // The ScatterNode's "scatter variable mapping" links to the innerGraph's scatter variable input Node:
        validatedOuterGraph.scatterNode.scatterVariableInnerGraphInputNodes.head eq validatedInnerGraph.x_scatterCollectionInput should be(true)

        // The ScatterNode's output port links to the inner graph's GraphOutputNode:
        validatedOuterGraph.scatterNode.outputMapping.toList match {
          case (port @ ScatterGathererPort(womType, outputToGather, _)) :: Nil =>
            port.name should be("foo.out")
            womType should be(WomArrayType(WomStringType))
            outputToGather eq validatedInnerGraph.foo_out_innerOutput should be(true)
          case other => fail("Expected exactly one output to be gathered in this scatter but got:" + other.mkString("\n", "\n", "\n"))
        }
      }

      val outer = validateOuterGraph
      val inner = validateInnerGraph(outer)
      validateConnections(outer, inner)
    }
  }

  it should "convert inputs into inner graph inputs" in {
    val scatterTest =
      """task foo {
        |  Int i
        |  command { ... }
        |  output {
        |    String str_out = i
        |    Int int_out = i + 1
        |  }
        |}
        |
        |workflow scatter_test {
        |  Int x = 5
        |  Int y = 6
        |  Int z = 7
        |  scatter (s in [x, y]) {
        |    call foo { input: i = z }
        |  }
        |}""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(scatterTest, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    val scatterTestGraph = namespace.workflow.toWomWorkflowDefinition.map(_.graph)

    scatterTestGraph match {
      case Valid(g) => validateGraph(g)
      case Invalid(errors) => fail(s"Unable to build wom version of scatter foo from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateGraph(workflowGraph: Graph) = {

      // Three inputs, a scatter node, the expression node for the scatter collection, and two outputs:
      workflowGraph.nodes.size should be(7)

      // Find that scatter:
      workflowGraph.nodes.collectFirst {
        case s: ScatterNode => s
      }.getOrElse(fail("Resulting graph did not contain a ScatterNode"))


    }
  }

  it should "convert unsatisfied call inputs in scatters into outer graph inputs" in {
    val scatterTest =
      """task foo {
        |  Int i
        |  Int j
        |  command { ... }
        |  output {
        |    String str_out = i
        |    Int int_out = i + 1
        |  }
        |}
        |
        |workflow scatter_test {
        |  Int x = 5
        |  scatter (s in range(x)) {
        |    call foo { input: i = s } # nb: foo.j is missing!
        |  }
        |}
        |""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(scatterTest, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    val scatterTestGraph = namespace.workflow.toWomWorkflowDefinition.map(_.graph)

    scatterTestGraph match {
      case Valid(g) => validateGraph(g)
      case Invalid(errors) => fail(s"Unable to build wom version of scatter foo from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateGraph(workflowGraph: Graph) = {

      // Find the inputs:
      val inputNodes: Set[ExternalGraphInputNode] = workflowGraph.nodes.filterByType[RequiredGraphInputNode]
      inputNodes.map {_.localName} should be(Set("foo.j"))
      inputNodes.map {_.identifier.fullyQualifiedName.value} should be(Set("scatter_test.foo.j"))

      // Find that scatter:
      val scatterNode = workflowGraph.nodes.collectFirst {
        case s: ScatterNode => s
      }.getOrElse(fail("Resulting graph did not contain a ScatterNode"))

      val scatterInnerInputs: Set[ExternalGraphInputNode] = scatterNode.innerGraph.nodes.filterByType[ExternalGraphInputNode]
      scatterInnerInputs map {_.identifier.fullyQualifiedName.value} should be(Set("scatter_test.foo.j"))
      val scatterInnerItemInput: Set[OuterGraphInputNode] = scatterNode.innerGraph.nodes.filterByType[OuterGraphInputNode]
      scatterInnerItemInput map {_.localName} should be(Set("s"))

      // Find the outputs:
      val outputNodes = workflowGraph.nodes.collect {
        case output: GraphOutputNode => output
      }
      outputNodes map { on => (on.localName, on.womType) } should be(Set(("foo.int_out", WomArrayType(WomIntegerType)), ("foo.str_out", WomArrayType(WomStringType))))

    }
  }

}
