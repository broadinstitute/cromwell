package wdl.transforms.biscayne.ast2wdlom

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements.{ExpressionElement, _}
import wdl.transforms.biscayne.ast2wdlom.WdlFileToWdlomSpec._
import wom.types._
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements._
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
          graphElements = Set(CallElement("in_n_out", None, None, Some(CallBodyElement(Vector(KvPair("total", IdentifierLookup("n")), KvPair("amount", IdentifierLookup("more"))))))),
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
      ),
    "afters" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "afters",
          inputsSection = None,
          graphElements = Set(
            CallElement("foo", None, None, Some(CallBodyElement(Vector(KvPair("i", ExpressionElement.PrimitiveLiteralExpressionElement(WomInteger(5))))))),
            CallElement("foo", Some("foo2"), Some("foo"), Some(CallBodyElement(Vector(KvPair("i", ExpressionElement.PrimitiveLiteralExpressionElement(WomInteger(6)))))))
          ),
          outputsSection = None,
          metaSection = None,
          parameterMetaSection = None)),
        tasks = Vector(TaskDefinitionElement(
          name = "foo",
          inputsSection = Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", None)))),
          declarations = Vector.empty,
          outputsSection = None,
          commandSection = CommandSectionElement(Vector(CommandSectionLine(Vector(
            StringCommandPartElement("cat \"hello "),
            PlaceholderCommandPartElement(IdentifierLookup("i"), PlaceholderAttributeSet.empty),
            StringCommandPartElement("\" > /tmp/helloFile")
          )))),
          runtimeSection = None,
          metaSection = None,
          parameterMetaSection = None
        ))
      ),
    "biscayne_escaping" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "escapes",
          inputsSection = None,
          graphElements = Set(
            IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "backslash", StringExpression(Vector(StringLiteral(" "), BackslashEscape, StringLiteral(" ")))),
            IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "dash", StringExpression(Vector(StringLiteral(" "), SingleQuoteEscape, StringLiteral(" ")))),
            IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "quote", StringExpression(Vector(StringLiteral(" "), DoubleQuoteEscape, StringLiteral(" ")))),
            IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "n", StringExpression(Vector(StringLiteral(" "), NewlineEscape, StringLiteral(" ")))),
            IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "t", StringExpression(Vector(StringLiteral(" "), TabEscape, StringLiteral(" ")))),
            IntermediateValueDeclarationElement(
              PrimitiveTypeElement(WomStringType),
              "octal_hello",
              StringExpression(Vector(UnicodeCharacterEscape(104), UnicodeCharacterEscape(101), UnicodeCharacterEscape(108), UnicodeCharacterEscape(108), UnicodeCharacterEscape(111)))
            ),
            IntermediateValueDeclarationElement(
              PrimitiveTypeElement(WomStringType),
              "hex_hello",
              StringExpression(Vector(UnicodeCharacterEscape(104), UnicodeCharacterEscape(101), UnicodeCharacterEscape(108), UnicodeCharacterEscape(108), UnicodeCharacterEscape(111)))
            ),
            IntermediateValueDeclarationElement(
              PrimitiveTypeElement(WomStringType),
              "unicode_hello",
              StringExpression(Vector(
                UnicodeCharacterEscape(104),
                UnicodeCharacterEscape(101),
                UnicodeCharacterEscape(108),
                UnicodeCharacterEscape(108),
                UnicodeCharacterEscape(111)
              ))
            )
          ),
          outputsSection = None,
          metaSection = None,
          parameterMetaSection = None
        )),
        tasks = Vector.empty)
  )
}
