package wdl.draft3.transforms.wdlom2wom

import cats.instances.either._
import better.files.File
import common.collections.EnhancedCollections._
import common.transforms.CheckedAtoB
import org.scalatest.{Assertion, FlatSpec, Matchers, Succeeded}
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.ast2wdlom._
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.CommandPartElement.StringCommandPartElement
import wdl.model.draft3.elements.ExpressionElement.StringLiteral
import wom.callable.Callable.{FixedInputDefinition, OptionalInputDefinition}
import wom.callable.MetaValueElement.{MetaValueElementBoolean, MetaValueElementObject}
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
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
      val converter: CheckedAtoB[File, WomBundle] = fileToAst andThen astToFileElement.map(fe => FileElementToWomBundleInputs(fe, "{}", List.empty, List.empty)) andThen fileElementToWomBundle

      converter.run(testCase) match {
        case Right(bundle) => validators(testName).apply(bundle)
        case Left(errors) =>
          val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WOM bundle: $formattedErrors")
      }
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
    "input_values" -> anyWomWillDo,
    "gap_in_command" -> anyWomWillDo,
    "nio_file" -> validateNioFile
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
        taskA.inputs.filter(_.isInstanceOf[FixedInputDefinition]).map(_.name).toSet should be(Set("rld", "__world1", "__world2"))
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
}
