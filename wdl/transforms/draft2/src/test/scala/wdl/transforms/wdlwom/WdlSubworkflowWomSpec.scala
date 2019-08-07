package wdl.transforms.wdlwom

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model.{Draft2ImportResolver, WdlNamespace, WdlNamespaceWithWorkflow}
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import wom.types.{WomArrayType, WomIntegerType, WomMaybeEmptyArrayType, WomStringType}
import wdl.transforms.draft2.wdlom2wom._
import wom.ResolvedImportRecord

class WdlSubworkflowWomSpec extends FlatSpec with Matchers {

  behavior of "WdlNamespaces with subworkflows"

  it should "support WDL to WOM conversion of subworkflow calls" in {
    val outerWdl =
      """import "import_me.wdl" as import_me
        |
        |workflow outer {
        |  Int x
        |  call import_me.inner as inner { input: i = x }
        |  output {
        |    Array[String] out = [inner.out]
        |  }
        |}""".stripMargin

    val innerWdl =
      """task foo {
        |  Int i
        |  command {
        |    echo ${i}
        |  }
        |  output {
        |    String out = read_string(stdout())
        |  }
        |}
        |
        |workflow inner {
        |  Int i
        |  call foo { input: i = i }
        |  output {
        |    String out = foo.out
        |    Int x = 5050
        |  }
        |}
      """.stripMargin


    def innerResolver: Draft2ImportResolver = str => Draft2ResolvedImportBundle(innerWdl, ResolvedImportRecord(str))

    val namespace = WdlNamespace.loadUsingSource(
      workflowSource = outerWdl,
      resource = None,
      importResolver = Some(Seq(innerResolver))).get.asInstanceOf[WdlNamespaceWithWorkflow]

    val outerWorkflowGraph = namespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).map(_.graph)

    outerWorkflowGraph match {
      case Valid(g) => validateOuter(g)
      case Invalid(errors) => fail(s"Unable to build wom version of workflow with subworkflow from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateOuter(workflowGraph: Graph) = {
      // One input, x
      workflowGraph.inputNodes.map(_.localName) should be(Set("x"))
      val calls = workflowGraph.calls
      calls.map(_.localName) should be(Set("inner"))

      // One workflow call, "inner"
      val innerCall = calls.head.asInstanceOf[WorkflowCallNode]
      innerCall.localName should be("inner")
      innerCall.identifier.fullyQualifiedName.value should be("outer.inner")
      innerCall.upstream.head.asInstanceOf[ExpressionNode].inputPorts.map(_.upstream.graphNode) should be(Set(workflowGraph.inputNodes.head))

      // One output, "out"
      workflowGraph.outputNodes.map(_.localName) should be(Set("out"))
      workflowGraph.outputNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("outer.out"))
      workflowGraph.outputNodes.head.womType should be(WomMaybeEmptyArrayType(WomStringType))
      workflowGraph.outputNodes.foreach(_.upstream should be(Set(innerCall)))

      validateInner(innerCall.callable.innerGraph)
    }

    def validateInner(innerGraph: Graph) = {
      innerGraph.inputNodes.map(_.localName) should be(Set("i"))
      val calls = innerGraph.calls
      calls.map(_.localName) should be(Set("foo"))

      val fooCall = calls.head.asInstanceOf[CommandCallNode]
      fooCall.upstream.head.asInstanceOf[ExpressionNode].inputPorts.map(_.upstream.graphNode) should be(Set(innerGraph.inputNodes.head))

      innerGraph.outputNodes.map(_.localName) should be(Set("out", "x"))
      innerGraph.outputNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("inner.out", "inner.x"))
      val outOutput = innerGraph.outputNodes.find(_.localName == "out").get
      val xOutput = innerGraph.outputNodes.find(_.localName == "x").get

      outOutput.womType should be(WomStringType)
      outOutput.upstream should be(Set(fooCall))

      xOutput.womType should be(WomIntegerType)
      xOutput.upstream should be(Set.empty)
    }
  }

  it should "support WDL to WOM conversion of subworkflows in scatters" in {
    val outerWdl =
      """import "import_me.wdl" as import_me
        |
        |workflow outer {
        |  Array[Int] xs
        |  scatter (x in xs) {
        |    call import_me.inner as inner { input: i = x }
        |  }
        |  output {
        |    Array[String] outs = inner.out
        |  }
        |}""".stripMargin

    val innerWdl =
      """task foo {
        |  Int i
        |  command {
        |    echo ${i}
        |  }
        |  output {
        |    String out = read_string(stdout())
        |  }
        |}
        |
        |workflow inner {
        |  Int i
        |  call foo { input: i = i }
        |  output {
        |    String out = foo.out
        |  }
        |}
      """.stripMargin


    def innerResolver: Draft2ImportResolver = str => Draft2ResolvedImportBundle(innerWdl, ResolvedImportRecord(str))

    val namespace = WdlNamespace.loadUsingSource(
      workflowSource = outerWdl,
      resource = None,
      importResolver = Some(Seq(innerResolver))).get.asInstanceOf[WdlNamespaceWithWorkflow]

    val outerWorkflowGraph = namespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).map(_.graph)

    outerWorkflowGraph match {
      case Valid(g) => validateOuter(g)
      case Invalid(errors) => fail(s"Unable to build wom version of workflow with subworkflow from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateOuter(workflowGraph: Graph) = {
      workflowGraph.inputNodes.map(_.localName) should be(Set("xs"))

      val scatter = workflowGraph.scatters.head
      val scatterCollectionNode = workflowGraph.nodes.collectFirst({ case e: ExpressionNode if e.localName == "x" => e }).get
      scatter.upstream should be(Set(scatterCollectionNode))
      scatter.outputPorts.map(_.name) should be(Set("inner.out"))
      scatter.outputPorts.head.womType should be(WomArrayType(WomStringType))

      workflowGraph.outputNodes.map(_.localName) should be(Set("outs"))
      workflowGraph.outputNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("outer.outs"))
      workflowGraph.outputNodes.head.womType should be(WomArrayType(WomStringType))
      workflowGraph.outputNodes.foreach(_.upstream should be(Set(scatter)))
    }
  }

  it should "support conversion of sub workflows with identical names" in {
    val outerWdl =
      """import "import_me.wdl" as import_me
        |
        |workflow twin {
        |  Int i = 5
        |  call import_me.twin { input: i = i }
        |  output {
        |    String outs = twin.out
        |  }
        |}""".stripMargin

    val innerWdl =
      """task foo {
        |  Int i
        |  command {
        |    echo ${i}
        |  }
        |  output {
        |    String out = read_string(stdout())
        |  }
        |}
        |
        |workflow twin {
        |  Int i
        |  call foo { input: i = i }
        |  output {
        |    String out = foo.out
        |  }
        |}
      """.stripMargin

    def innerResolver: Draft2ImportResolver = str => Draft2ResolvedImportBundle(innerWdl, ResolvedImportRecord(str))

    val namespace = WdlNamespace.loadUsingSource(
      workflowSource = outerWdl,
      resource = None,
      importResolver = Some(Seq(innerResolver))).get.asInstanceOf[WdlNamespaceWithWorkflow]

    val outerWorkflowGraph = namespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).map(_.graph)

    outerWorkflowGraph match {
      case Valid(g) => validateOuter(g)
      case Invalid(errors) => fail(s"Unable to build wom version of workflow with subworkflow from WDL: ${errors.toList.mkString("\n", "\n", "\n")}")
    }

    def validateOuter(workflowGraph: Graph) = {
      // One input, x
      workflowGraph.inputNodes.size should be(1)
      workflowGraph.inputNodes foreach {
        case opt: OptionalGraphInputNodeWithDefault =>
          opt.localName should be("i")
          opt.womType should be(WomIntegerType)
        case other => fail("Unexpected input node to outer graph: " + other.localName)
      }
      val calls = workflowGraph.calls
      calls.map(_.localName) should be(Set("twin"))

      // One workflow call, "twin"
      val innerCall = calls.head.asInstanceOf[WorkflowCallNode]
      innerCall.localName should be("twin")
      innerCall.identifier.fullyQualifiedName.value should be("twin.twin")

      // One output, "out"
      workflowGraph.outputNodes.map(_.localName) should be(Set("outs"))
      workflowGraph.outputNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("twin.outs"))
      workflowGraph.outputNodes.head.womType should be(WomStringType)
      workflowGraph.outputNodes.foreach(_.upstream should be(Set(innerCall)))

      validateInner(innerCall.callable.innerGraph)
    }

    def validateInner(innerGraph: Graph) = {
      innerGraph.inputNodes.map(_.localName) should be(Set("i"))
      val calls = innerGraph.calls
      calls.map(_.localName) should be(Set("foo"))

      val fooCall = calls.head.asInstanceOf[CommandCallNode]
      fooCall.upstream.head.asInstanceOf[ExpressionNode].inputPorts.map(_.upstream.graphNode) should be(Set(innerGraph.inputNodes.head))

      innerGraph.outputNodes.map(_.localName) should be(Set("out"))
      innerGraph.outputNodes.map(_.identifier.fullyQualifiedName.value) should be(Set("twin.out"))
      val outOutput = innerGraph.outputNodes.find(_.localName == "out").get

      outOutput.womType should be(WomStringType)
      outOutput.upstream should be(Set(fooCall))
    }
  }

}
