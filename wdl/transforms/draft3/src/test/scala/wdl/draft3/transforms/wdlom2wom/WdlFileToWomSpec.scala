package wdl.draft3.transforms.wdlom2wom

import cats.instances.either._
import better.files.File
import common.collections.EnhancedCollections._
import common.transforms.CheckedAtoB
import org.scalatest.{Assertion, FlatSpec, Matchers, Succeeded}
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.ast2wdlom._
import wdl.transforms.base.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.CommandPartElement.StringCommandPartElement
import wdl.model.draft3.elements.ExpressionElement.StringLiteral
import wdl.transforms.base.wdlom2wom._
import wom.callable.Callable.{FixedInputDefinitionWithDefault, OptionalInputDefinition}
import wom.callable.MetaValueElement._
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wom.graph.expression.{ExposedExpressionNode, TaskCallInputExpressionNode}
import wom.graph.{ScatterNode, WorkflowCallNode}
import wom.types._

class WdlFileToWomSpec extends FlatSpec with Matchers {
  behavior of "WDL File to WOM"

  val testCases = File("wdl/transforms/draft3/src/test/cases")

  it should "be set up for testing" in {
    testCases.exists shouldBe true
    testCases.list.nonEmpty shouldBe true
  }

  testCases.list.filter(x => x.isRegularFile && x.extension.contains(".wdl")) foreach { testCase =>

    val fileName = testCase.name
    val testName = testCase.name.split("\\.").head

    val itShouldString = s"create a valid WOM object for $fileName"
    val testOrIgnore: (=>Any) => Unit = if (testCase.name.endsWith(".ignored.wdl") || testCase.name.endsWith(".nowom.wdl")) {
      (it should itShouldString).ignore _
    } else {
      (it should itShouldString).in _
    }

    testOrIgnore {
      val converter: CheckedAtoB[File, WomBundle] = fileToAst andThen wrapAst andThen astToFileElement.map(fe => FileElementToWomBundleInputs(fe, "{}", convertNestedScatterToSubworkflow = true, List.empty, List.empty, workflowDefinitionElementToWomWorkflowDefinition, taskDefinitionElementToWomTaskDefinition)) andThen fileElementToWomBundle

      converter.run(testCase) match {
        case Right(bundle) => validators(testName).apply(bundle)
        case Left(errors) =>
          val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WOM bundle: $formattedErrors")
      }
    }
  }

  // There is a scatter within a scatter
  //
  // scatter(a in indices) {
  //   scatter(b in indices) {
  //     Int x = a + b
  //   }
  // }
  //
  it should "be able to leave nested scatters intact" in {
    val converter: CheckedAtoB[File, WomBundle] = fileToAst andThen wrapAst andThen astToFileElement.map(fe => FileElementToWomBundleInputs(fe, "{}", convertNestedScatterToSubworkflow = false, List.empty, List.empty, workflowDefinitionElementToWomWorkflowDefinition, taskDefinitionElementToWomTaskDefinition)) andThen fileElementToWomBundle

    val twoLevelScatterFile = File("wdl/transforms/draft3/src/test/cases/two_level_scatter.wdl")

    converter.run(twoLevelScatterFile) match {
      case Right(bundle) =>
        val wf = bundle.primaryCallable.get.asInstanceOf[WorkflowDefinition]
        val graph = wf.innerGraph

        // get the top scatter node
        graph.scatters.size shouldBe(1)
        val topScatter : ScatterNode = graph.scatters.toVector.head
        val wfCalls = graph.allNodes.filterByType[WorkflowCallNode]

        // don't generate any sub-workflows
        wfCalls.size shouldBe(0)

        // there should be one scatter inside the top scatter
        val innerGraph = topScatter.innerGraph
        innerGraph.scatters.size shouldBe(1)
        Succeeded

      case Left(errors) =>
        val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
        fail(s"Failed to create WOM bundle: $formattedErrors")
    }
  }


  it should "split a nested scatter into a toplevel scatter, and a bottom sub-workflow" in {
    val converter: CheckedAtoB[File, WomBundle] = fileToAst andThen wrapAst andThen astToFileElement.map(fe => FileElementToWomBundleInputs(fe, "{}", convertNestedScatterToSubworkflow = true, List.empty, List.empty, workflowDefinitionElementToWomWorkflowDefinition, taskDefinitionElementToWomTaskDefinition)) andThen fileElementToWomBundle

    val twoLevelScatterFile = File("wdl/transforms/draft3/src/test/cases/two_level_scatter.wdl")

    converter.run(twoLevelScatterFile) match {
      case Right(bundle) =>
        val wf = bundle.primaryCallable.get.asInstanceOf[WorkflowDefinition]
        val graph = wf.innerGraph

        // There should be just one scatter.
        graph.scatters.size shouldBe(1)
        val wfCalls = graph.allNodes.filterByType[WorkflowCallNode]

        // There should be a call to a generated sub-workflow in the graph
        wfCalls.size shouldBe(1)
        Succeeded
      case Left(errors) =>
        val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
        fail(s"Failed to create WOM bundle: $formattedErrors")
    }
  }


  private val validators: Map[String, WomBundle => Assertion] = Map(
    "declaration_chain" -> anyWomWillDo,
    "empty_workflow" -> anyWomWillDo,
    "input_expressions" -> anyWomWillDo,
    "input_types" -> anyWomWillDo,
    "input_values" -> anyWomWillDo,
    "passthrough_workflow" -> anyWomWillDo,
    "scatter_var_member_access" -> anyWomWillDo,
    "nested_conditionals" -> anyWomWillDo,
    "simple_first_test" -> anyWomWillDo,
    "static_value_workflow" -> anyWomWillDo,
    "nested_struct" -> anyWomWillDo,
    "struct_definition" -> validateStructDefinitionWom,
    "simple_scatter" -> anyWomWillDo,
    "ogin_scatter" -> anyWomWillDo,
    "nested_scatter" -> anyWomWillDo,
    "two_level_scatter" -> anyWomWillDo,
    "simple_conditional" -> anyWomWillDo,
    "lots_of_nesting" -> anyWomWillDo,
    "standalone_task" -> anyWomWillDo,
    "simple_task" -> validateTaskDefinitionWom,
    "lots_of_nesting" -> anyWomWillDo,
    "taskless_engine_functions" -> anyWomWillDo,
    "no_input_no_output_workflow" -> anyWomWillDo,
    "command_syntaxes" -> validateCommandSyntaxes,
    "standalone_task" -> anyWomWillDo,
    "task_with_metas" -> anyWomWillDo,
    "task_with_metas2" -> validateMetaSection,
    "input_values" -> anyWomWillDo,
    "gap_in_command" -> anyWomWillDo,
    "nio_file" -> validateNioFile,
    "same_named_inputs_priority" -> validateSameNamedInputsPriority,
    "cmd_whitespace_none" -> anyWomWillDo,
    "cmd_strip_common_spaces" -> anyWomWillDo,
    "cmd_whitespace_tabs" -> anyWomWillDo,
    "cmd_strip_common_tabs" -> anyWomWillDo,
    "cmd_whitespace_spaces" -> anyWomWillDo,
    "string_escaping" -> anyWomWillDo
  )

  private def anyWomWillDo(b: WomBundle): Assertion = Succeeded

  private def validateStructDefinitionWom(b: WomBundle): Assertion = {
    val wfDef: WorkflowDefinition = (b.allCallables.values.toSet.filterByType[WorkflowDefinition]: Set[WorkflowDefinition]).head
    b.typeAliases.keySet shouldBe Set("FooStruct")
    val structOutputType = (wfDef.graph.outputNodes.map(_.womType).filterByType[WomCompositeType]: Set[WomCompositeType]).head

    structOutputType.typeMap shouldBe Map(
      "simple" -> WomIntegerType,
      "complex" -> WomPairType(WomArrayType(WomIntegerType), WomMapType(WomStringType, WomBooleanType))
    )
  }

  private def validateTaskDefinitionWom(b: WomBundle): Assertion = {
    val taskDef: CallableTaskDefinition = (b.allCallables.values.toSet.filterByType[CallableTaskDefinition]: Set[CallableTaskDefinition]).head
    taskDef.name shouldBe "simple"
    taskDef.commandTemplate(Map.empty) shouldBe List(WdlomWomStringCommandPart(StringCommandPartElement("echo Hello World ")))
  }

  private def validateCommandSyntaxes(b: WomBundle): Assertion = {
    b.allCallables.size should be(2)
    b.allCallables.get("a")match {
      case Some(taskA) =>
        taskA.inputs.filter(_.isInstanceOf[FixedInputDefinitionWithDefault]).map(_.name).toSet should be(Set("rld", "__world1", "__world2"))
        taskA.inputs.filter(_.isInstanceOf[OptionalInputDefinition]).map(_.name).toSet should be(Set("world1", "world2"))
        taskA.inputs.map(_.name).toSet should be(Set("rld", "__world1", "__world2", "world1", "world2"))
        taskA.outputs.map(_.name).toSet should be(Set("out"))
        taskA.asInstanceOf[CallableTaskDefinition].runtimeAttributes.attributes("docker").asInstanceOf[WdlomWomExpression].expressionElement should be(StringLiteral("ubuntu:latest"))
      case None => fail("Expected a task called 'a'")
    }
    b.allCallables.get("b") match {
      case Some(taskB) =>
        taskB.inputs.map(_.name) should be(Seq("world"))
        taskB.outputs.map(_.name) should be(Seq("out"))
        taskB.asInstanceOf[CallableTaskDefinition].runtimeAttributes.attributes("docker").asInstanceOf[WdlomWomExpression].expressionElement should be(StringLiteral("ubuntu:latest"))
      case None => fail("Expected a task called 'b'")
    }
  }

  private def validateNioFile(b: WomBundle): Assertion = {
    b.allCallables.size should be(1)
    b.allCallables.get("nio_file") match {
      case None => fail("No callable found 'nio_file'")
      case Some(nioFileTask) =>
        // Plain old input:
        nioFileTask.inputs.find(_.name == "f").get.parameterMeta should be(Some(MetaValueElementObject(Map("localization_optional" -> MetaValueElementBoolean(true)))))
        // Input based on upstream:
        nioFileTask.inputs.find(_.name == "g").get.parameterMeta should be(Some(MetaValueElementObject(Map("localization_optional" -> MetaValueElementBoolean(true)))))
        // Optional input:
        nioFileTask.inputs.find(_.name == "h").get.parameterMeta should be(Some(MetaValueElementObject(Map("localization_optional" -> MetaValueElementBoolean(true)))))
    }
  }

  private def validateSameNamedInputsPriority(b: WomBundle): Assertion = {
    val nodes = b.primaryCallable.get.asInstanceOf[WorkflowDefinition].graph.nodes

    val exposedExpressionNode = nodes.filter(_.isInstanceOf[ExposedExpressionNode]).head
    val callInputs = nodes.filter(_.isInstanceOf[TaskCallInputExpressionNode]).toList

    callInputs.size shouldBe 4

    // Usually, we use node names to check whether in-memory graphs are linked together correctly.
    // Obviously, in this regression test for a same-named-node bug that won't work - instead,
    // we use identity (the hash code) to differentiate multiple nodes with the same name.
    callInputs(0).inputPorts.head.upstream should be theSameInstanceAs exposedExpressionNode.outputPorts.head
    callInputs(1).inputPorts.head.upstream should be theSameInstanceAs exposedExpressionNode.outputPorts.head
    callInputs(2).inputPorts.head.upstream should be theSameInstanceAs exposedExpressionNode.outputPorts.head
    callInputs(3).inputPorts.head.upstream should be theSameInstanceAs exposedExpressionNode.outputPorts.head
  }

  private def validateMetaSection(b: WomBundle): Assertion = {
    val task = b.primaryCallable.get.asInstanceOf[CallableTaskDefinition]

    task.meta should be (Map("author" -> MetaValueElementString("John Doe"),
                             "email" -> MetaValueElementString("john.doe@yahoo.com"),
                             "b" -> MetaValueElementBoolean(true),
                             "zipcode" -> MetaValueElementInteger(94043),
                             "f" -> MetaValueElementFloat(1.3),
                             "numbers" -> MetaValueElementArray(Vector(MetaValueElementInteger(1),
                                                                       MetaValueElementInteger(2),
                                                                       MetaValueElementInteger(3))),
                             "extras" -> MetaValueElementObject(Map("house" -> MetaValueElementString("With porch"),
                                                                    "cat" -> MetaValueElementString("Lucy")))
                         ))
  }
}
