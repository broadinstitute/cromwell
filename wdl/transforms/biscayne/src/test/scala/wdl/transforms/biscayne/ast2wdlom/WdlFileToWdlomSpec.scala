package wdl.transforms.biscayne.ast2wdlom

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements._
import wdl.transforms.biscayne.ast2wdlom.WdlFileToWdlomSpec._
import wom.types._
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement._

import wom.values.WomInteger

class WdlFileToWdlomSpec extends FlatSpec with Matchers {

  behavior of "WDL File to WDLOM"

  val testCases = File("wdl/transforms/biscayne/src/test/cases")

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
    "simple_first_test" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "order",
          inputsSection = Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "n", Some(PrimitiveLiteralExpressionElement(WomInteger(4)))),
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "more", Some(StringLiteral("more")))))),
          graphElements = Set(CallElement("in_n_out", None, Some(CallBodyElement(Vector(KvPair("total", IdentifierLookup("n")), KvPair("amount", IdentifierLookup("more"))))))),
          outputsSection = Some(OutputsSectionElement(Vector(
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "out", IdentifierMemberAccess("in_n_out", "out", List.empty))))),
          metaSection = None,
          parameterMetaSection = None)),
        tasks = Vector(TaskDefinitionElement(
          name = "in_n_out",
          inputsSection = Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "total", None),
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "amount", None)))),
          declarations = Vector.empty,
          outputsSection = Some(OutputsSectionElement(Vector(
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "out", Add(ReadInt(StdoutElement), PrimitiveLiteralExpressionElement(WomInteger(1))))))),
          commandSection = CommandSectionElement(Vector(CommandSectionLine(Vector(
            StringCommandPartElement("echo "),
            PlaceholderCommandPartElement(IdentifierLookup("total"), PlaceholderAttributeSet.empty),
            StringCommandPartElement(" ")
          )))),
          runtimeSection = None,
          metaSection = None,
          parameterMetaSection = None
        ))
      )
  )
}
