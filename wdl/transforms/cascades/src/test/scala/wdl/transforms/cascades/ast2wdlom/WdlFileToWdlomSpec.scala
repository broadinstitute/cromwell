package wdl.transforms.cascades.ast2wdlom

import better.files.File
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.elements._
import wdl.transforms.cascades.ast2wdlom.WdlFileToWdlomSpec._
import wom.SourceFileLocation
import wom.types._
import wom.values.{WomBoolean, WomInteger}

class WdlFileToWdlomSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "WDL File to WDLOM"

  val testCases = File("wdl/transforms/cascades/src/test/cases")

  it should "be set up for testing" in {
    testCases.exists shouldBe true
    testCases.list.nonEmpty shouldBe true
  }

  testCases.list.filter(x => x.isRegularFile && x.extension.contains(".wdl")) foreach { testCase =>
    val fileName = testCase.name
    val testName = testCase.name.split("\\.").head

    val itShouldString = s"create the correct Element structure for $fileName"
    val testOrIgnore: (=> Any) => Unit = if (fileName.endsWith(".ignored.wdl")) {
      (it should itShouldString).ignore _
    } else {
      (it should itShouldString).in _
    }

    testOrIgnore {

      val expected = expectations.getOrElse(testName, fail(s"No Element expectation defined for $testName"))
      fileToFileElement.run(testCase) match {
        case Right(actual) => actual shouldBe expected
        case Left(errors) =>
          val formattedErrors =
            errors.toList.mkString(System.lineSeparator(), System.lineSeparator(), System.lineSeparator())
          fail(s"Failed to create WDLOM:$formattedErrors")
      }

    }
  }
}

object WdlFileToWdlomSpec {

  val expectations: Map[String, FileElement] = Map(
    "no_input_no_output_workflow" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(
          WorkflowDefinitionElement(
            name = "no_input_no_output",
            inputsSection = None,
            graphElements = Set(
              CallElement("no_inputs", Some("noi1"), Vector.empty, None, Some(SourceFileLocation(7))),
              CallElement("no_inputs", Some("noi5"), Vector.empty, None, Some(SourceFileLocation(21))),
              CallElement("no_inputs", Some("noi2"), Vector.empty, None, Some(SourceFileLocation(11))),
              CallElement("no_inputs", Some("noi3"), Vector.empty, None, Some(SourceFileLocation(13))),
              CallElement("no_inputs", Some("noi4"), Vector.empty, None, Some(SourceFileLocation(17))),
              CallElement("no_inputs", None, Vector.empty, None, Some(SourceFileLocation(9)))
            ),
            outputsSection = None,
            metaSection = None,
            parameterMetaSection = None,
            sourceLocation = Some(SourceFileLocation(6))
          )
        ),
        tasks = Vector(
          TaskDefinitionElement(
            name = "no_inputs",
            inputsSection = None,
            declarations = Vector.empty,
            outputsSection = None,
            commandSection =
              CommandSectionElement(List(CommandSectionLine(Vector(StringCommandPartElement("echo Hello World "))))),
            runtimeSection = None,
            metaSection = None,
            parameterMetaSection = None,
            sourceLocation = Some(SourceFileLocation(27))
          )
        )
      ),
    "simple_first_test" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(
          WorkflowDefinitionElement(
            name = "order",
            inputsSection = Some(
              InputsSectionElement(
                Vector(
                  InputDeclarationElement(PrimitiveTypeElement(WomIntegerType),
                                          "n",
                                          Some(PrimitiveLiteralExpressionElement(WomInteger(4)))
                  ),
                  InputDeclarationElement(PrimitiveTypeElement(WomStringType), "more", Some(StringLiteral("more")))
                )
              )
            ),
            graphElements = Set(
              CallElement(
                "in_n_out",
                None,
                Vector.empty,
                Some(
                  CallBodyElement(
                    Vector(KvPair("total", IdentifierLookup("n")), KvPair("amount", IdentifierLookup("more")))
                  )
                ),
                Some(SourceFileLocation(19))
              )
            ),
            outputsSection = Some(
              OutputsSectionElement(
                Vector(
                  OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType),
                                           "out",
                                           IdentifierMemberAccess("in_n_out", "out", List.empty)
                  )
                )
              )
            ),
            metaSection = None,
            parameterMetaSection = None,
            sourceLocation = Some(SourceFileLocation(14))
          )
        ),
        tasks = Vector(
          TaskDefinitionElement(
            name = "in_n_out",
            inputsSection = Some(
              InputsSectionElement(
                Vector(InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "total", None),
                       InputDeclarationElement(PrimitiveTypeElement(WomStringType), "amount", None)
                )
              )
            ),
            declarations = Vector.empty,
            outputsSection = Some(
              OutputsSectionElement(
                Vector(
                  OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType),
                                           "out",
                                           Add(ReadInt(StdoutElement), PrimitiveLiteralExpressionElement(WomInteger(1)))
                  )
                )
              )
            ),
            commandSection = CommandSectionElement(
              Vector(
                CommandSectionLine(
                  Vector(
                    StringCommandPartElement("echo "),
                    PlaceholderCommandPartElement(IdentifierLookup("total"), PlaceholderAttributeSet.empty),
                    StringCommandPartElement(" ")
                  )
                )
              )
            ),
            runtimeSection = None,
            metaSection = None,
            parameterMetaSection = None,
            sourceLocation = Some(SourceFileLocation(3))
          )
        )
      ),
    "afters" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(
          WorkflowDefinitionElement(
            name = "afters",
            inputsSection = None,
            graphElements = Set(
              CallElement(
                "foo",
                None,
                Vector.empty,
                Some(
                  CallBodyElement(
                    Vector(KvPair("i", ExpressionElement.PrimitiveLiteralExpressionElement(WomInteger(5))))
                  )
                ),
                Some(SourceFileLocation(4))
              ),
              CallElement(
                "foo",
                Some("foo2"),
                Vector("foo"),
                Some(
                  CallBodyElement(
                    Vector(KvPair("i", ExpressionElement.PrimitiveLiteralExpressionElement(WomInteger(6))))
                  )
                ),
                Some(SourceFileLocation(5))
              )
            ),
            outputsSection = None,
            metaSection = None,
            parameterMetaSection = None,
            sourceLocation = Some(SourceFileLocation(3))
          )
        ),
        tasks = Vector(
          TaskDefinitionElement(
            name = "foo",
            inputsSection = Some(
              InputsSectionElement(Vector(InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", None)))
            ),
            declarations = Vector.empty,
            outputsSection = None,
            commandSection = CommandSectionElement(
              Vector(
                CommandSectionLine(
                  Vector(
                    StringCommandPartElement("cat \"hello "),
                    PlaceholderCommandPartElement(IdentifierLookup("i"), PlaceholderAttributeSet.empty),
                    StringCommandPartElement("\" > /tmp/helloFile")
                  )
                )
              )
            ),
            runtimeSection = None,
            metaSection = None,
            parameterMetaSection = None,
            sourceLocation = Some(SourceFileLocation(8))
          )
        )
      ),
    "cascades_escaping" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(
          WorkflowDefinitionElement(
            name = "escapes",
            inputsSection = None,
            graphElements = Set(
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "backslash",
                StringExpression(Vector(StringLiteral(" "), BackslashEscape, StringLiteral(" ")))
              ),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "n",
                StringExpression(Vector(StringLiteral(" "), NewlineEscape, StringLiteral(" ")))
              ),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "t",
                StringExpression(Vector(StringLiteral(" "), TabEscape, StringLiteral(" ")))
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType),
                                                  "q1",
                                                  StringExpression(
                                                    Vector(StringLiteral("leading text "),
                                                           DoubleQuoteEscape,
                                                           StringLiteral(" trailing text")
                                                    )
                                                  )
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "q2", StringLiteral("\"")),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "q3",
                StringExpression(Vector(StringLiteral("  "), DoubleQuoteEscape, StringLiteral("  ")))
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType),
                                                  "q4",
                                                  StringExpression(
                                                    Vector(StringLiteral("leading text "),
                                                           SingleQuoteEscape,
                                                           StringLiteral(" trailing text")
                                                    )
                                                  )
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "q5", StringLiteral("'")),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "q6",
                StringExpression(Vector(StringLiteral("  "), SingleQuoteEscape, StringLiteral("  ")))
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType),
                                                  "sq1",
                                                  StringExpression(
                                                    Vector(StringLiteral("leading text "),
                                                           DoubleQuoteEscape,
                                                           StringLiteral(" trailing text")
                                                    )
                                                  )
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "sq2", StringLiteral("\"")),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "sq3",
                StringExpression(Vector(StringLiteral("  "), DoubleQuoteEscape, StringLiteral("  ")))
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType),
                                                  "sq4",
                                                  StringExpression(
                                                    Vector(StringLiteral("leading text "),
                                                           SingleQuoteEscape,
                                                           StringLiteral(" trailing text")
                                                    )
                                                  )
              ),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "sq5", StringLiteral("'")),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "sq6",
                StringExpression(Vector(StringLiteral("  "), SingleQuoteEscape, StringLiteral("  ")))
              ),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "octal_hello",
                StringExpression(
                  Vector(UnicodeCharacterEscape(104),
                         UnicodeCharacterEscape(101),
                         UnicodeCharacterEscape(108),
                         UnicodeCharacterEscape(108),
                         UnicodeCharacterEscape(111)
                  )
                )
              ),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "hex_hello",
                StringExpression(
                  Vector(UnicodeCharacterEscape(104),
                         UnicodeCharacterEscape(101),
                         UnicodeCharacterEscape(108),
                         UnicodeCharacterEscape(108),
                         UnicodeCharacterEscape(111)
                  )
                )
              ),
              IntermediateValueDeclarationElement(
                PrimitiveTypeElement(WomStringType),
                "unicode_hello",
                StringExpression(
                  Vector(
                    UnicodeCharacterEscape(104),
                    UnicodeCharacterEscape(101),
                    UnicodeCharacterEscape(108),
                    UnicodeCharacterEscape(108),
                    UnicodeCharacterEscape(111)
                  )
                )
              )
            ),
            outputsSection = None,
            metaSection = None,
            parameterMetaSection = None,
            sourceLocation = Some(SourceFileLocation(3))
          )
        ),
        tasks = Vector.empty
      ),
    "struct_literal" -> FileElement(
      imports = Vector(),
      structs = Vector(
        StructElement("Plant",
                      Vector(StructEntryElement("color", PrimitiveTypeElement(WomStringType)),
                             StructEntryElement("tasty", PrimitiveTypeElement(WomBooleanType))
                      )
        ),
        StructElement(
          "Animal",
          Vector(StructEntryElement("name", PrimitiveTypeElement(WomStringType)),
                 StructEntryElement("isGood", OptionalTypeElement(PrimitiveTypeElement(WomBooleanType)))
          )
        )
      ),
      workflows = Vector(
        WorkflowDefinitionElement(
          "struct_literal",
          None,
          Set(
            CallElement(
              "test_struct_parsing",
              None,
              Vector(),
              Some(
                CallBodyElement(
                  Vector(
                    KvPair("standard_plant_input",
                           StructLiteral("Plant",
                                         Map("color" -> StringLiteral("green"),
                                             "tasty" -> PrimitiveLiteralExpressionElement(WomBoolean(true))
                                         )
                           )
                    ),
                    KvPair("standard_animal_input",
                           StructLiteral("Animal",
                                         Map("name" -> StringLiteral("mittens"),
                                             "isGood" -> PrimitiveLiteralExpressionElement(WomBoolean(false))
                                         )
                           )
                    )
                  )
                )
              ),
              Some(SourceFileLocation(35))
            )
          ),
          None,
          None,
          None,
          Some(SourceFileLocation(34))
        )
      ),
      tasks = Vector(
        TaskDefinitionElement(
          "test_struct_parsing",
          Some(
            InputsSectionElement(
              Vector(
                InputDeclarationElement(TypeAliasElement("Plant"), "standard_plant_input", None),
                InputDeclarationElement(TypeAliasElement("Animal"), "standard_animal_input", None)
              )
            )
          ),
          Vector(),
          Some(
            OutputsSectionElement(
              Vector(
                OutputDeclarationElement(TypeAliasElement("Plant"),
                                         "standard_plant_forwarded",
                                         IdentifierLookup("standard_plant_input")
                ),
                OutputDeclarationElement(TypeAliasElement("Animal"),
                                         "standard_animal_forwarded",
                                         IdentifierLookup("standard_animal_input")
                ),
                OutputDeclarationElement(
                  TypeAliasElement("Plant"),
                  "plant_output_literal",
                  StructLiteral("Plant",
                                Map("color" -> StringLiteral("red"),
                                    "tasty" -> PrimitiveLiteralExpressionElement(WomBoolean(true))
                                )
                  )
                )
              )
            )
          ),
          CommandSectionElement(
            List(CommandSectionLine(Vector(StringCommandPartElement("""echo "all dogs are good""""))))
          ),
          Some(RuntimeAttributesSectionElement(Vector(KvPair("docker", StringLiteral("ubuntu:latest"))))),
          None,
          None,
          Some(SourceFileLocation(13))
        )
      )
    )
  )
}
