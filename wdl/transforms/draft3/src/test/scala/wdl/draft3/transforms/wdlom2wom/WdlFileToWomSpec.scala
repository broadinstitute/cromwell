package wdl.draft3.transforms.wdlom2wom

import cats.instances.either._
import better.files.File
import common.collections.EnhancedCollections._
import common.transforms.CheckedAtoB
import org.scalatest.{Assertion, FlatSpec, Matchers, Succeeded}
import wdl.draft3.transforms.parsing._
import wdl.draft3.transforms.ast2wdlom._
import wom.callable.WorkflowDefinition
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
      val converter: CheckedAtoB[File, WomBundle] = fileToAst andThen astToFileElement.map(fe => FileElementAndImportResolvers(fe, List.empty)) andThen fileElementToWomBundle

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
    "simpleFirstTest" -> anyWomWillDo,
    "static_value_workflow" -> anyWomWillDo,
    "nested_struct" -> anyWomWillDo,
    "struct_definition" -> validateStructDefinitionWom

  )

  private def anyWomWillDo(b: WomBundle) = Succeeded

  private def validateStructDefinitionWom(b: WomBundle): Assertion = {
    val wfDef: WorkflowDefinition = (b.callables.filterByType[WorkflowDefinition]: Set[WorkflowDefinition]).head
    b.typeAliases.keySet shouldBe Set("FooStruct")
    val structOutputType = (wfDef.graph.outputNodes.map(_.womType).filterByType[WomCompositeType]: Set[WomCompositeType]).head

    structOutputType.typeMap shouldBe Map(
      "simple" -> WomIntegerType,
      "complex" -> WomPairType(WomArrayType(WomIntegerType), WomMapType(WomStringType, WomBooleanType))
    )
  }
}
