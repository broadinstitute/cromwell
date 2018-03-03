package wdl.draft3.transforms.ast2wdlom

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements._
import wdl.draft3.transforms.ast2wdlom.WdlFileToWdlomSpec._
import wom.types._
import wdl.draft3.transforms.ast2wdlom.ExpressionSet._
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement._
import wom.values.{WomBoolean, WomInteger}

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
            PairTypeElement(ArrayTypeElement(PrimitiveTypeElement(WomIntegerType)),MapTypeElement(PrimitiveTypeElement(WomStringType),PrimitiveTypeElement(WomBooleanType)))
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
            "complex" -> PairLiteral(ArrayLiteral(Vector(PrimitiveLiteralExpressionElement(WomInteger(5)))),MapLiteral(Map(StringLiteral("t") -> PrimitiveLiteralExpressionElement(WomBoolean(true)))))
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
        imports = List.empty,
        structs = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
         name = "order",
         inputsSection = Some(InputsSectionElement(Vector(
           InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "n", Some(PrimitiveLiteralExpressionElement(WomInteger(4)))),
           InputDeclarationElement(PrimitiveTypeElement(WomStringType), "more", Some(StringLiteral("more")))))),
         graphElements = Set(CallElement("in_n_out", None, Some(CallBodyElement(Vector(KvPair("total", IdentifierLookup("n")), KvPair("amount", IdentifierLookup("more"))))))),
         outputsSection = None,
          metaSection = None,
          parameterMetaSection = None)),
        tasks = Vector(TaskDefinitionElement(
          name = "in_n_out",
          inputsSection = Some(InputsSectionElement(Vector(
            InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "total", None),
            InputDeclarationElement(PrimitiveTypeElement(WomStringType), "amount", None)))),
          outputsSection = Some(OutputsSectionElement(Vector(
            OutputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "out", Add(ReadString(StdoutElement), PrimitiveLiteralExpressionElement(WomInteger(1))))))),
          commandSection = CommandSectionElement(Vector(StringCommandPartElement(" echo "), PlaceholderCommandPartElement(IdentifierLookup("total")), StringCommandPartElement(" "))),
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
        workflows = Vector.empty,
        tasks = Vector(
          TaskDefinitionElement(
            name = "standalone",
            inputsSection = Some(InputsSectionElement(Vector(InputDeclarationElement(PrimitiveTypeElement(WomStringType), "bar", None)))),
            outputsSection = Some(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomStringType), "out", IdentifierLookup("bar"))))),
            commandSection = CommandSectionElement(Vector(StringCommandPartElement("\n    echo "), PlaceholderCommandPartElement(IdentifierLookup("bar")), StringCommandPartElement("\n  "))),
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
            outputsSection = Some(OutputsSectionElement(Vector.empty)),
            commandSection = CommandSectionElement(Vector(StringCommandPartElement(" echo Hello World "))),
            runtimeSection = None,
            metaSection = Some(MetaSectionElement(Vector(
                                               (MetaKvPair("author", MetaValueElement.MString("John Doe"))),
                                               (MetaKvPair("email", MetaValueElement.MString("john.doe@yahoo.com")))
                                             ))),
            parameterMetaSection = Some(ParameterMetaSectionElement(
                                          Vector(
                                            (MetaKvPair("a", MetaValueElement.MString("just an integer"))),
                                            (MetaKvPair("b", MetaValueElement.MString("an important parameter")))
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
            outputsSection = None,
            commandSection = CommandSectionElement(Vector(StringCommandPartElement(" echo Hello World "))),
            runtimeSection = None,
            metaSection = None,
            parameterMetaSection = None
          )
        )
      )
  )
}
