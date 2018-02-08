package wdl.transforms.wdlwom

import cats.data.Validated.{Invalid, Valid}
import common.collections.EnhancedCollections._
import org.scalatest.{FlatSpec, Matchers}
import wom.graph._
import wom.transforms.WomWorkflowDefinitionMaker
import wdl.WdlWorkflow

class WdlAliasWomSpec(implicit workflowDefinitionMaker: WomWorkflowDefinitionMaker[WdlWorkflow]) extends FlatSpec with Matchers {

  behavior of "WdlNamespaces with aliased calls"

  it should "get successfully converted into WOM" in {
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
        |  call foo as foo1
        |  call foo as foo2
        |}""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(conditionalTest, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
    val conditionalTestGraph = namespace.workflow.toWomWorkflowDefinition.map(_.graph)

    conditionalTestGraph match {
      case Valid(g) => validateGraph(g)
      case Invalid(errors) => fail(s"Unable to build wom version of conditional foo from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateGraph(workflowGraph: Graph) = {

      val inputNodes: Set[ExternalGraphInputNode] = workflowGraph.nodes.filterByType[ExternalGraphInputNode]
      inputNodes.map(_.localName) should be(Set("foo1.i", "foo2.i"))
      inputNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("conditional_test.foo1.i", "conditional_test.foo2.i"))

      val callNodes: Set[CallNode] = workflowGraph.nodes.filterByType[CallNode]
      callNodes.map(_.localName) should be(Set("foo1", "foo2"))
      callNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("conditional_test.foo1", "conditional_test.foo2"))

      val outputNodes: Set[GraphOutputNode] = workflowGraph.nodes.filterByType[GraphOutputNode]
      outputNodes.map(_.localName) should be(Set("foo1.out", "foo2.out"))
      outputNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("conditional_test.foo1.out", "conditional_test.foo2.out"))

      workflowGraph.nodes.size should be(6)
    }
  }
}
