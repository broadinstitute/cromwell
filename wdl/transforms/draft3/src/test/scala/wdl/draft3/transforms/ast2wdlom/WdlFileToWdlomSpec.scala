package wdl.draft3.transforms.ast2wdlom

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements._
import wdl.draft3.transforms.ast2wdlom.WdlFileToWdlomSpec._
import wom.types._
import wdl.draft3.transforms.ast2wdlom.ExpressionSet._
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement._
import wom.values.{WomBoolean, WomFloat, WomInteger}
import wdl.model.draft3.elements.MetaValueElement._

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
        structs = Vector.empty,
        workflows = List(WorkflowDefinitionElement(
          name = "empty",
          inputsSection = None,
          graphElements = Set.empty,
          outputsSection = None,
          metaSection = None,
          parameterMetaSection = None)),
        tasks = List.empty),
    "struct_definition" -> FileElement(
      imports = Vector(),
      structs = Vector(StructElement(
        name = "FooStruct",
        entries = Vector(
          StructEntryElement(
            "simple",
            PrimitiveTypeElement(WomIntegerType)),
          StructEntryElement(
            "complex",
            PairTypeElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), MapTypeElement(PrimitiveTypeElement(WomStringType), PrimitiveTypeElement(WomBooleanType)))
          )
        )
      )),
      workflows = Vector(WorkflowDefinitionElement(
        name = "struct_definition",
        inputsSection = None,
        graphElements = Set(),
        outputsSection = Some(OutputsSectionElement(Vector(OutputDeclarationElement(
          typeElement = TypeAliasElement("FooStruct"),
          name = "myFoo",
          expression = ObjectLiteral(Map(
            "simple" -> PrimitiveLiteralExpressionElement(WomInteger(5)),
            "complex" -> PairLiteral(ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(5)))), MapLiteral(Map(StringLiteral("t") -> PrimitiveLiteralExpressionElement(WomBoolean(true)))))
          ))
        )))),
        metaSection = None,
        parameterMetaSection = None)),
      tasks = Vector.empty
    ),
    "input_types" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "input_types",
          inputsSection = Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", None),
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "s", None),
            InputDeclarationElement(PrimitiveTypeElement(WomFloatType), "float", None),
            InputDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b", None),
            InputDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "file", None),
            InputDeclarationElement(ObjectTypeElement, "o", None),
            InputDeclarationElement(OptionalTypeElement(PrimitiveTypeElement(WomIntegerType)), "maybe_i", None),
            InputDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomStringType)), "array_s", None),
            InputDeclarationElement(MapTypeElement(PrimitiveTypeElement(WomIntegerType), PrimitiveTypeElement(WomStringType)), "map_is", None),
            InputDeclarationElement(
              ArrayTypeElement(
                OptionalTypeElement(
                  PairTypeElement(PrimitiveTypeElement(WomStringType), PrimitiveTypeElement(WomIntegerType)))),
              "lotsa_nesting_array", None)
          ))),
          graphElements = Set.empty,
          outputsSection = None,
          metaSection = None,
          parameterMetaSection = None)),
        tasks = Vector.empty),
    "input_values" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "input_values",
          inputsSection = Some(InputsSectionElement(
            inputDeclarations = Vector(
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", Some(intLiteral)),
              InputDeclarationElement(PrimitiveTypeElement(WomStringType), "s", Some(stringLiteral)),
              InputDeclarationElement(PrimitiveTypeElement(WomStringType), "placeholder", Some(stringPlaceholderExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomStringType), "placeholder2", Some(stringPlaceholderExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomFloatType), "f", Some(floatLiteral)),
              InputDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b", Some(booleanLiteral))
            )
          )),
          graphElements = Set.empty,
          outputsSection = None,
          metaSection = None,
          parameterMetaSection = None)),
        tasks = List.empty),
    "input_expressions" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "input_expressions",
          inputsSection = Some(InputsSectionElement(
            inputDeclarations = Vector(
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "ten", Some(addExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "zero", Some(subtractExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "twentyfive", Some(multiplyExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "one", Some(divideExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "zeroagain", Some(remainderExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "tenagain", Some(tenVariableLookup)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "pair_expression_member_access", Some(pairExpressionMemberAccess)),
              InputDeclarationElement(PrimitiveTypeElement(WomBooleanType), "unary_expressions", Some(unaryExpressions)),
              InputDeclarationElement(PrimitiveTypeElement(WomBooleanType), "comparisons", Some(comparisonExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "variableLookupMemberAccesses", Some(chainIdentifierAccess)),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "expressionMemberAccesses", Some(chainPairAccess)),
              InputDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "is", Some(arrayOfIs)),
              InputDeclarationElement(ObjectTypeElement, "object_literal", Some(objectLiteralExpression)),
              InputDeclarationElement(
                MapTypeElement(PrimitiveTypeElement(WomStringType), PrimitiveTypeElement(WomIntegerType)),
                "map_literal",
                Some(mapLiteralExpression)
              ),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "ternaryIf", Some(ternaryIfExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomStringType), "string_read", Some(ReadString(StdoutElement))),
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "zipped", Some(zippedExpression)),
              InputDeclarationElement(PrimitiveTypeElement(WomStringType), "subbed", Some(subbedExpression))
            )
          )),
          graphElements = Set.empty,
          outputsSection = None,
          metaSection = None,
          parameterMetaSection = None)),
        tasks = Vector.empty),
    "passthrough_workflow" ->
      FileElement(
        imports = Vector(),
        structs = Vector.empty,
        workflows = Vector(
          WorkflowDefinitionElement(
            name = "foo",
            inputsSection = Some(InputsSectionElement(Vector(
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "x", None)
            ))),
            graphElements = Set(IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "y", IdentifierLookup("x"))),
            outputsSection = Some(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "z", IdentifierLookup("y"))))),
            metaSection = None,
            parameterMetaSection = None)),
        tasks = Vector()),
    "declaration_chain" ->
      FileElement(
        imports = Vector(),
        structs = Vector.empty,
        workflows = Vector(
          WorkflowDefinitionElement(
            name = "foo",
            inputsSection = None,
            graphElements = Set(
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "y", IdentifierLookup("x")),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "a", intLiteral),
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "x", IdentifierLookup("a"))
            ),
            outputsSection = Some(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "z", IdentifierLookup("y"))))),
            metaSection = None,
            parameterMetaSection = None)),
        tasks = Vector()),
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
        ))),
    "static_value_workflow" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "foo",
          inputsSection = None,
          graphElements = Set.empty,
          outputsSection = Some(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "y", PrimitiveLiteralExpressionElement(WomInteger(3)))))),
          metaSection = None,
          parameterMetaSection = None)),
        tasks = Vector.empty),
    "standalone_task" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement("standalone_task", None, Set.empty, None, None, None)),
        tasks = Vector(
          TaskDefinitionElement(
            name = "standalone",
            inputsSection = Some(InputsSectionElement(Vector(InputDeclarationElement(PrimitiveTypeElement(WomStringType), "bar", None)))),
            declarations = Vector.empty,
            outputsSection = Some(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomStringType), "out", IdentifierLookup("bar"))))),
            commandSection = CommandSectionElement(Vector(
              CommandSectionLine(Vector(
                StringCommandPartElement("echo "),
                PlaceholderCommandPartElement(IdentifierLookup("bar"), PlaceholderAttributeSet.empty)
              ))
            )),
            runtimeSection = Some(RuntimeAttributesSectionElement(Vector(KvPair("docker", StringLiteral("someFakeDockerRuntime"))))),
            metaSection = None,
            parameterMetaSection = None))
      ),
    "task_with_metas" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector.empty,
        tasks = Vector(
          TaskDefinitionElement(
            name = "task_with_metas",
            inputsSection = Some(InputsSectionElement(Vector.empty)),
            declarations = Vector.empty,
            outputsSection = Some(OutputsSectionElement(Vector.empty)),
            commandSection = CommandSectionElement(Vector(CommandSectionLine(Vector(
              StringCommandPartElement("echo Hello World ")
            )))),
            runtimeSection = None,
            metaSection = Some(MetaSectionElement(
              Map("author" -> MetaValueElementString("John Doe"),
                "email" -> MetaValueElementString("john.doe@yahoo.com"))
            )),
            parameterMetaSection = Some(ParameterMetaSectionElement(
              Map("a" -> MetaValueElementString("just an integer"),
                "b" -> MetaValueElementString("an important parameter"),
                "x" -> MetaValueElementArray(Vector(MetaValueElementString("A"),
                  MetaValueElementString("B"),
                  MetaValueElementString("C"))),
                "y" -> MetaValueElementArray(Vector(MetaValueElementInteger(1),
                  MetaValueElementInteger(2),
                  MetaValueElementInteger(3))),
                "yf" -> MetaValueElementArray(Vector(MetaValueElementFloat(1.1),
                  MetaValueElementFloat(2.9),
                  MetaValueElementFloat(3.14))),
                "z" -> MetaValueElementObject(Map("k1" -> MetaValueElementInteger(1),
                  "k2" -> MetaValueElementInteger(2),
                  "k3" -> MetaValueElementInteger(3)))
              )))
          ))),
    "no_input_no_output_workflow" ->
      FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "no_input_no_output",
          inputsSection = None,
          graphElements = Set(CallElement("no_inputs", None, None)),
          outputsSection = None,
          metaSection = None,
          parameterMetaSection = None)),
        tasks = Vector(
          TaskDefinitionElement(
            name = "no_inputs",
            inputsSection = None,
            declarations = Vector.empty,
            outputsSection = None,
            commandSection = CommandSectionElement(Vector(CommandSectionLine(Vector(StringCommandPartElement("echo Hello World "))))),
            runtimeSection = None,
            metaSection = None,
            parameterMetaSection = None
          )
        )
      ),
    "nested_struct" ->
      FileElement(
        imports = Vector(),
        structs = Vector(
          StructElement(
            name = "A",
            entries = Vector(
              StructEntryElement("i", PrimitiveTypeElement(WomIntegerType)),
              StructEntryElement("f", PrimitiveTypeElement(WomFloatType)))
          ),
          StructElement(
            name = "B",
            entries = Vector(
              StructEntryElement("a", TypeAliasElement("A")),
              StructEntryElement("i", PrimitiveTypeElement(WomIntegerType)),
              StructEntryElement("f", PrimitiveTypeElement(WomFloatType)))
          )),
        workflows = Vector(WorkflowDefinitionElement(
          name = "nested_struct",
          inputsSection = None,
          graphElements = Set(
            IntermediateValueDeclarationElement(
              TypeAliasElement("B"),
              "b",
              ObjectLiteral(Map(
                "a" -> ObjectLiteral(Map(
                  "i" -> PrimitiveLiteralExpressionElement(WomInteger(5)),
                  "f" -> PrimitiveLiteralExpressionElement(WomFloat(5.5)))),
                "i" -> PrimitiveLiteralExpressionElement(WomInteger(6)),
                "f" -> PrimitiveLiteralExpressionElement(WomFloat(6.6)))))),
          outputsSection = Some(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomFloatType), "f", IdentifierMemberAccess("b", "a", Vector("f")))))),
          metaSection = None,
          parameterMetaSection = None
        )),
        tasks = Vector()
      ),
    "simple_scatter" -> FileElement(
      imports = Vector.empty,
      structs = Vector.empty,
      workflows = Vector(WorkflowDefinitionElement(
        name = "simple_scatter",
        inputsSection = None,
        graphElements = Set(
          IntermediateValueDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "indices", ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(1)), PrimitiveLiteralExpressionElement(WomInteger(2)), PrimitiveLiteralExpressionElement(WomInteger(3))))),
          ScatterElement(
            scatterName = "ScatterAt6_11",
            scatterExpression = IdentifierLookup("indices"),
            scatterVariableName = "i",
            graphElements = Vector(
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "j", Add(IdentifierLookup("i"), PrimitiveLiteralExpressionElement(WomInteger(10))))
            )
          )
        ),
        outputsSection = Some(
          OutputsSectionElement(Vector(OutputDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "js", IdentifierLookup("j"))))
        ),
        metaSection = None,
        parameterMetaSection = None
      )),
      tasks = Vector.empty
    ),
    "ogin_scatter" -> FileElement(
      imports = Vector.empty,
      structs = Vector.empty,
      workflows = Vector(WorkflowDefinitionElement(
        name = "ogin_scatter",
        inputsSection = None,
        graphElements = Set(
          IntermediateValueDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "indices", ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(1)), PrimitiveLiteralExpressionElement(WomInteger(2)), PrimitiveLiteralExpressionElement(WomInteger(3))))),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "ogin_me", PrimitiveLiteralExpressionElement(WomInteger(10))),
          ScatterElement(
            scatterName = "ScatterAt8_11",
            scatterExpression = IdentifierLookup("indices"),
            scatterVariableName = "i",
            graphElements = Vector(
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "j", Add(IdentifierLookup("i"), IdentifierLookup("ogin_me")))
            )
          )
        ),
        outputsSection = Some(
          OutputsSectionElement(Vector(OutputDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "js", IdentifierLookup("j"))))
        ),
        metaSection = None,
        parameterMetaSection = None
      )),
      tasks = Vector.empty
    ),
    "nested_scatter" -> FileElement(
      Vector(),
      Vector(),
      Vector(WorkflowDefinitionElement(
        "nested_scatter",
        None,
        Set(
          IntermediateValueDeclarationElement(
            ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)),
            "indices",
            ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(1)), PrimitiveLiteralExpressionElement(WomInteger(2)), PrimitiveLiteralExpressionElement(WomInteger(3))))
          ),
          IntermediateValueDeclarationElement(
            PrimitiveTypeElement(WomIntegerType),
            "y",
            PrimitiveLiteralExpressionElement(WomInteger(55))
          ),
          ScatterElement(
            scatterName = "ScatterAt8_11",
            IdentifierLookup("indices"),
            "a",
            Vector(
              ScatterElement(
                scatterName = "ScatterAt9_13",
                IdentifierLookup("indices"), "b",
                Vector(
                  IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "x", Add(IdentifierLookup("a"), IdentifierLookup("b"))),
                  ScatterElement(
                    scatterName = "ScatterAt11_15",
                    IdentifierLookup("indices"),
                    "c",
                    Vector(
                      IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "j", Add(Add(Add(IdentifierLookup("a"), IdentifierLookup("b")), IdentifierLookup("c")), IdentifierLookup("x"))))),
                  ScatterElement(
                    scatterName = "ScatterAt14_15",
                    IdentifierLookup("j"),
                    "d",
                    Vector(IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "k", Add(IdentifierLookup("d"), IdentifierLookup("y"))))
                  )
                )
              )
            )
          )
        ),
        Some(
          OutputsSectionElement(Vector(OutputDeclarationElement(ArrayTypeElement(ArrayTypeElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)))), "ks", IdentifierLookup("k"))))
        ),
        None,
        None)
      ),
      Vector()
    ),
    "simple_conditional" -> FileElement(
      imports = Vector.empty,
      structs = Vector.empty,
      workflows = Vector(WorkflowDefinitionElement(
        name = "simple_conditional",
        inputsSection = None,
        graphElements = Set(
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomBooleanType), "bool", PrimitiveLiteralExpressionElement(WomBoolean(true))),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", PrimitiveLiteralExpressionElement(WomInteger(5))),
          IfElement(
            conditionExpression = IdentifierLookup("bool"),
            graphElements = Vector(
              IntermediateValueDeclarationElement(PrimitiveTypeElement(WomIntegerType), "j", Add(IdentifierLookup("i"), PrimitiveLiteralExpressionElement(WomInteger(10))))
            )
          )
        ),
        outputsSection = Some(
          OutputsSectionElement(Vector(OutputDeclarationElement(OptionalTypeElement(PrimitiveTypeElement(WomIntegerType)), "j_maybe", IdentifierLookup("j"))))
        ),
        metaSection = None,
        parameterMetaSection = None
      )),
      tasks = Vector.empty
    ),
    "lots_of_nesting" -> FileElement(
      imports = Vector.empty,
      structs = Vector.empty,
      workflows = Vector(WorkflowDefinitionElement(
        name = "lots_of_nesting",
        inputsSection = None,
        graphElements = Set(
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b0", PrimitiveLiteralExpressionElement(WomBoolean(true))),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b1", PrimitiveLiteralExpressionElement(WomBoolean(true))),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b2", PrimitiveLiteralExpressionElement(WomBoolean(true))),
          IntermediateValueDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "i0s", Range(PrimitiveLiteralExpressionElement(WomInteger(2)))),
          IntermediateValueDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "i1s", Range(PrimitiveLiteralExpressionElement(WomInteger(2)))),
          IntermediateValueDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "i2s", Range(PrimitiveLiteralExpressionElement(WomInteger(2)))),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "s0", StringLiteral("hello")),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "s1", StringLiteral("world")),
          IfElement(
            conditionExpression = IdentifierLookup("b0"),
            graphElements = Vector(
              ScatterElement(
                scatterName = "ScatterAt17_13",
                scatterExpression = IdentifierLookup("i0s"),
                scatterVariableName = "i0",
                graphElements = Vector(
                  IfElement(
                    conditionExpression = IdentifierLookup("b1"),
                    graphElements = Vector(
                      ScatterElement(
                        scatterName = "ScatterAt19_17",
                        scatterExpression = IdentifierLookup("i1s"),
                        scatterVariableName = "i1",
                        graphElements = Vector(
                          IfElement(
                            conditionExpression = IdentifierLookup("b2"),
                            graphElements = Vector(
                              ScatterElement(
                                "ScatterAt21_21",
                                IdentifierLookup("i2s"),
                                "i2",
                                Vector(
                                  IntermediateValueDeclarationElement(
                                    PrimitiveTypeElement(WomStringType),
                                    "s",
                                    Add(IdentifierLookup("s0"), IdentifierLookup("s1"))
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        ),
        outputsSection = Some(OutputsSectionElement(Vector(
          OutputDeclarationElement(
            typeElement = OptionalTypeElement(ArrayTypeElement(OptionalTypeElement(ArrayTypeElement(OptionalTypeElement(ArrayTypeElement(PrimitiveTypeElement(WomStringType))))))),
            name = "s_out",
            expression = IdentifierLookup("s"))))),
        metaSection = None,
        parameterMetaSection = None
      )),
      tasks = Vector.empty
    ),
    "simple_task" -> FileElement(
        imports = Vector.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement("simple_task", None, Set.empty, None, None, None)),
        tasks = Vector(
          TaskDefinitionElement(
            name = "simple",
            inputsSection = None,
            declarations = Vector.empty,
            outputsSection = None,
            commandSection = CommandSectionElement(Vector(CommandSectionLine(Vector(StringCommandPartElement("echo Hello World "))))),
            runtimeSection = None,
            metaSection = None,
            parameterMetaSection = None))
    ),
    "default_input_overrides" -> null,
    "taskless_engine_functions" -> FileElement(
      imports = Vector.empty,
      structs = Vector.empty,
      workflows = Vector(WorkflowDefinitionElement(
        name = "taskless_engine_functions",
        inputsSection = None,
        graphElements = Set(
          IntermediateValueDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "ints", ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(1)), PrimitiveLiteralExpressionElement(WomInteger(2))))),
          IntermediateValueDeclarationElement(NonEmptyTypeElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType))), "definitelyInts", ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(1)), PrimitiveLiteralExpressionElement(WomInteger(2))))),
          IntermediateValueDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomStringType)), "strings", ArrayLiteral(Vector(StringLiteral("a"), StringLiteral("b")))),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "filepath", StringLiteral("gs://not/a/real/file.txt")),
          IntermediateValueDeclarationElement(ArrayTypeElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType))), "matrix",
            ArrayLiteral(Vector(
              ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(1)), PrimitiveLiteralExpressionElement(WomInteger(0)))),
              ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(1)), PrimitiveLiteralExpressionElement(WomInteger(0))))
            ))
          ),
          IntermediateValueDeclarationElement(ArrayTypeElement(MapTypeElement(PrimitiveTypeElement(WomIntegerType), PrimitiveTypeElement(WomStringType))), "list_of_maps",
            ArrayLiteral(Vector(
              MapLiteral(Map(PrimitiveLiteralExpressionElement(WomInteger(1)) -> StringLiteral("one"), PrimitiveLiteralExpressionElement(WomInteger(2)) -> StringLiteral("two"))),
              MapLiteral(Map(PrimitiveLiteralExpressionElement(WomInteger(11)) -> StringLiteral("eleven"), PrimitiveLiteralExpressionElement(WomInteger(22)) -> StringLiteral("twenty-two")))
            ))
          ),
          IntermediateValueDeclarationElement(PrimitiveTypeElement(WomFloatType), "f", PrimitiveLiteralExpressionElement(WomFloat(1.024)))
        ),
        outputsSection = Some(
          OutputsSectionElement(Vector(
            OutputDeclarationElement(ArrayTypeElement(PairTypeElement(PrimitiveTypeElement(WomIntegerType),PrimitiveTypeElement(WomStringType))), "int_cross_string", Cross(IdentifierLookup("ints"),IdentifierLookup("strings"))),
            OutputDeclarationElement(ArrayTypeElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType))), "transposed_matrix", Transpose(IdentifierLookup("matrix"))),
            OutputDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "flattened_matrix", Flatten(IdentifierLookup("matrix"))),
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "matrix_length", Length(IdentifierLookup("matrix"))),
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "flattened_matrix_length", Length(IdentifierLookup("flattened_matrix"))),
            OutputDeclarationElement(ArrayTypeElement(PairTypeElement(PrimitiveTypeElement(WomIntegerType),PrimitiveTypeElement(WomStringType))), "flattened_map", Flatten(IdentifierLookup("list_of_maps"))),
            OutputDeclarationElement(PrimitiveTypeElement(WomStringType), "file_basename", Basename(IdentifierLookup("filepath"),None)),
            OutputDeclarationElement(PrimitiveTypeElement(WomStringType), "file_basename_extensionless", Basename(IdentifierLookup("filepath"),Some(StringLiteral(".txt")))),
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "f_floor",Floor(IdentifierLookup("f"))),
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "f_ceiling",Ceil(IdentifierLookup("f"))),
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "f_round",Round(IdentifierLookup("f"))),
            OutputDeclarationElement(
              ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)), "m1", IndexAccess(IdentifierLookup("matrix"),PrimitiveLiteralExpressionElement(WomInteger(1)))),
            OutputDeclarationElement(
              PrimitiveTypeElement(WomIntegerType),"m2", IndexAccess(
                IndexAccess(IdentifierLookup("matrix"),PrimitiveLiteralExpressionElement(WomInteger(1))),PrimitiveLiteralExpressionElement(WomInteger(1))))
          ))
        ),
        metaSection = None,
        parameterMetaSection = None
      )),
      tasks = Vector.empty
    ),
    "command_syntaxes" -> FileElement(
      imports = Vector.empty,
      structs = Vector.empty,
      workflows = Vector.empty,
      tasks = Vector(
        TaskDefinitionElement(
          name = "a",
          inputsSection = Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "world1", Some(StringExpression(Vector(StringLiteral("wo"), StringPlaceholder(IdentifierLookup("rld")))))),
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "world2", Some(StringExpression(Vector(StringLiteral("wo"), StringPlaceholder(IdentifierLookup("rld"))))))
          ))),
          declarations = Vector(
            IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "rld", StringLiteral("rld"))
          ),
          outputsSection = Some(OutputsSectionElement(Vector(
            OutputDeclarationElement(PrimitiveTypeElement(WomStringType), "out", ReadString(StdoutElement))
          ))),
          commandSection = CommandSectionElement(Vector(
            CommandSectionLine(Vector(
              StringCommandPartElement("echo "),
              PlaceholderCommandPartElement(StringLiteral("hello"), PlaceholderAttributeSet.empty),
              StringCommandPartElement(" "),
              PlaceholderCommandPartElement(IdentifierLookup("world1"), PlaceholderAttributeSet.empty)
            )),
            CommandSectionLine(Vector(
              StringCommandPartElement("echo goodbye "),
              PlaceholderCommandPartElement(IdentifierLookup("world2"), PlaceholderAttributeSet.empty)
            )),
            CommandSectionLine(Vector(
              StringCommandPartElement("echo "),
              PlaceholderCommandPartElement(IdentifierLookup("world1"),
              PlaceholderAttributeSet(Some("foo"), Some("--yes"), Some("--no"), Some(", ")))
            ))
          )),
          runtimeSection = Some(RuntimeAttributesSectionElement(Vector(
            KvPair("docker", StringLiteral("ubuntu:latest"))
          ))),
          metaSection = None,
          parameterMetaSection = None
        ),
        TaskDefinitionElement(
          name = "b",
          inputsSection = Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "world", Some(StringLiteral("world")))
          ))),
          declarations = Vector.empty,
          outputsSection = Some(OutputsSectionElement(Vector(
            OutputDeclarationElement(PrimitiveTypeElement(WomStringType), "out", ReadString(StdoutElement))
          ))),
          commandSection = CommandSectionElement(Vector(
            CommandSectionLine(Vector(
              StringCommandPartElement("echo hello ${world}")
            )),
            CommandSectionLine(Vector(
              StringCommandPartElement("echo goodbye "),
              PlaceholderCommandPartElement(IdentifierLookup("world"), PlaceholderAttributeSet.empty)
            )
            ))),
          runtimeSection = Some(RuntimeAttributesSectionElement(Vector(
            KvPair("docker", StringLiteral("ubuntu:latest"))
          ))),
          metaSection = None,
          parameterMetaSection = None
        )
      )
    ),
    "gap_in_command" -> FileElement(
      imports = Vector.empty,
      structs = Vector.empty,
      workflows = Vector(WorkflowDefinitionElement(
        "my_workflow",
        None,
        Set(CallElement("my_task",None,None)),
        None,
        None,
        None
      )),
      tasks = Vector(TaskDefinitionElement(
        "my_task",
        None,
        Vector(),
        Some(OutputsSectionElement(
          Vector(OutputDeclarationElement(ArrayTypeElement(PrimitiveTypeElement(WomStringType)),"lines",ReadLines(StdoutElement)))
        )),
        CommandSectionElement(Vector(CommandSectionLine(Vector(StringCommandPartElement("""    echo "hi""""))), CommandSectionLine(Vector()), CommandSectionLine(Vector(StringCommandPartElement("""    echo "bye""""))))),
        None,
        None,
        None))
    )
  )
}
