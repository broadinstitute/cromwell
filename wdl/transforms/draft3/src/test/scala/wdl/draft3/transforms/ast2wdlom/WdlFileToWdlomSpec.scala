package wdl.draft3.transforms.ast2wdlom

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements._
import wdl.draft3.transforms.ast2wdlom.WdlFileToWdlomSpec._
import wdl.model.draft3.elements.ExpressionElement._
import wom.types._
import wom.values.{WomBoolean, WomFloat, WomInteger, WomString}

class WdlFileToWdlomSpec extends FlatSpec with Matchers {

  behavior of "WDL File to WDLOM"

  val testCases = File("wdl/transforms/draft3/src/test/cases")

  it should "be set up for testing" in {
    testCases.exists shouldBe true
    testCases.list.nonEmpty shouldBe true
  }

  testCases.list.filter(x => x.isRegularFile && x.extension.contains(".wdl")) foreach { testCase =>

    val fileName = testCase.name
    val testName = testCase.name.split("\\.").head

    val itShouldString = s"create the correct Element structure for $fileName"
    val testOrIgnore: (=>Any) => Unit = if (fileName.endsWith(".ignored.wdl")) {
      (it should itShouldString).ignore _
    } else {
      (it should itShouldString).in _
    }

    testOrIgnore {

      val expected = expectations.getOrElse(testName, fail(s"No Element expectation defined for $testName"))
      fileToFileElement.run(testCase) match {
        case Right(actual) => actual shouldBe expected
        case Left(errors) =>
          val formattedErrors = errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WDLOM:$formattedErrors")
      }

    }
  }
}

object WdlFileToWdlomSpec {

  val expectations: Map[String, FileElement] = Map(
    "empty_workflow" ->
      FileElement(
        imports = List.empty,
        workflows = List(WorkflowDefinitionElement("empty", None, Vector.empty)),
        tasks = List.empty),
    "input_types" ->
      FileElement(
        imports = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          "input_types",
          Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", None),
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "s", None),
            InputDeclarationElement(PrimitiveTypeElement(WomFloatType), "f", None),
            InputDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b", None),
            InputDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "f", None),
            InputDeclarationElement(ObjectTypeElement, "o", None),
            InputDeclarationElement(OptionalTypeElement(PrimitiveTypeElement(WomIntegerType)), "maybe_i", None),
            InputDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomStringType)), "array_s", None),
            InputDeclarationElement(MapTypeElement(PrimitiveTypeElement(WomIntegerType), PrimitiveTypeElement(WomStringType)), "map_is", None),
            InputDeclarationElement(
              ArrayTypeElement(
                OptionalTypeElement(
                  PairTypeElement(PrimitiveTypeElement(WomStringType), PrimitiveTypeElement(WomIntegerType)))),
              "lotsa_nesting_array", None)
          ))), Vector.empty
        )),
        tasks = Vector.empty),
    "input_values" ->
      FileElement(
        imports = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "input_values",
          inputsSection = Some(InputsSectionElement(
            inputDeclarations = Vector(
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", Some(PrimitiveLiteralExpressionElement(WomInteger(5)))),
              InputDeclarationElement(PrimitiveTypeElement(WomStringType), "s", Some(PrimitiveLiteralExpressionElement(WomString("s")))),
              InputDeclarationElement(PrimitiveTypeElement(WomFloatType), "f", Some(PrimitiveLiteralExpressionElement(WomFloat(5.5)))),
              InputDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b", Some(PrimitiveLiteralExpressionElement(WomBoolean(true))))
            )
          )),
          outputsSection = Vector.empty)
        ),
        tasks = List.empty),
    "input_expressions" ->
      FileElement(
        imports = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "input_expressions",
          inputsSection = Some(InputsSectionElement(
            inputDeclarations = Vector(
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "four", Some(Add(
                left = PrimitiveLiteralExpressionElement(WomInteger(2)),
                right = PrimitiveLiteralExpressionElement(WomInteger(2)))))
            )
          )),
          outputsSection = Vector.empty)
        ),
        tasks = List.empty),
    "passthrough_workflow" ->
      FileElement(
        imports = List.empty,
        workflows = List.empty,
        tasks = List.empty),
    "simpleFirstTest" ->
      FileElement(
        imports = List.empty,
        workflows = List.empty,
        tasks = List.empty),
    "static_value_workflow" ->
      FileElement(
        imports = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement("foo", None, Vector(OutputsSectionElement(Vector(OutputElement("Int", "y", "3")))))),
        tasks = Vector.empty
      )
  )
}
