package wdl.draft3.transforms.ast2wdlom

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements._
import wdl.draft3.transforms.ast2wdlom.WdlFileToWdlomSpec._
import wom.types._
import wdl.draft3.transforms.ast2wdlom.ExpressionSet._
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement, ReadString, StdoutElement}
import wom.values.WomInteger

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
        workflows = List(WorkflowDefinitionElement("empty", None, None)),
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
          ))), None
        )),
        tasks = Vector.empty),
    "input_values" ->
      FileElement(
        imports = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement(
          name = "input_values",
          inputsSection = Some(InputsSectionElement(
            inputDeclarations = Vector(
              InputDeclarationElement(PrimitiveTypeElement(WomIntegerType), "i", Some(intLiteral)),
              InputDeclarationElement(PrimitiveTypeElement(WomStringType), "s", Some(stringLiteral)),
              InputDeclarationElement(PrimitiveTypeElement(WomFloatType), "f", Some(floatLiteral)),
              InputDeclarationElement(PrimitiveTypeElement(WomBooleanType), "b", Some(booleanLiteral))
            )
          )),
          outputsSection = None)
        ),
        tasks = List.empty),
    "input_expressions" ->
      FileElement(
        imports = Vector.empty,
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
          outputsSection = None)
        ),
        tasks = Vector.empty),
    "passthrough_workflow" ->
      FileElement(
        imports = Vector(),
        workflows = Vector(
          WorkflowDefinitionElement(
            name = "foo",
            Some(InputsSectionElement(Vector(InputDeclarationElement(PrimitiveTypeElement(WomIntegerType),"x",None)))),
            Some(OutputsSectionElement(Vector(DeclarationElement(PrimitiveTypeElement(WomIntegerType), "y", IdentifierLookup("x"))))))),
        tasks = Vector()),
    "simpleFirstTest" ->
      FileElement(
        imports = List.empty,
        workflows = List.empty,
        tasks = List.empty),
    "static_value_workflow" ->
      FileElement(
        imports = Vector.empty,
        workflows = Vector(WorkflowDefinitionElement("foo", None, Some(OutputsSectionElement(Vector(DeclarationElement(PrimitiveTypeElement(WomIntegerType), "y", PrimitiveLiteralExpressionElement(WomInteger(3)))))))),
        tasks = Vector.empty
      )
  )
}
