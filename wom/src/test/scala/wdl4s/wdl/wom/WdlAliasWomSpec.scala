package wdl.wom

import cats.data.Validated.{Invalid, Valid}
import lenthall.collections.EnhancedCollections._
import org.scalatest.{FlatSpec, Matchers}
import wdl.{WdlNamespace, WdlNamespaceWithWorkflow}
import wom.graph._

class WdlAliasWomSpec extends FlatSpec with Matchers {

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
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    val conditionalTestGraph = namespace.womExecutable.flatMap(_.graph)

    conditionalTestGraph match {
      case Valid(g) => validateGraph(g)
      case Invalid(errors) => fail(s"Unable to build wom version of conditional foo from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateGraph(workflowGraph: Graph) = {

      val inputNodes: Set[GraphInputNode] = workflowGraph.nodes.filterByType[GraphInputNode]
      inputNodes.map(_.name) should be(Set("conditional_test.foo1.i", "conditional_test.foo2.i"))

      val callNodes: Set[CallNode] = workflowGraph.nodes.filterByType[CallNode]
      callNodes.map(_.name) should be(Set("foo1", "foo2"))

      val outputNodes: Set[GraphOutputNode] = workflowGraph.nodes.filterByType[GraphOutputNode]
      outputNodes.map(_.name) should be(Set("foo1.out", "foo2.out"))

      workflowGraph.nodes.size should be(6)
    }
  }
}
