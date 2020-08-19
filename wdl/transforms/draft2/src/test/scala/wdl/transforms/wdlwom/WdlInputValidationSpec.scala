package wdl.transforms.wdlwom

import cats.data.NonEmptyList
import cats.syntax.either._
import common.Checked
import common.validation.Checked._
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor3}
import shapeless.Coproduct
import wdl.draft2.model.{WdlNamespace, WdlNamespaceWithWorkflow}
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wdl.transforms.draft2.wdlom2wom._
import wom.executable.Executable.ResolvedExecutableInputs
import wom.expression.NoIoFunctionSet
import wom.graph.Graph.ResolvedExecutableInput
import wom.graph.GraphNodePort.OutputPort
import wom.transforms.WomExecutableMaker.ops._
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import wom.types._
import wom.values._

class WdlInputValidationSpec extends AnyFlatSpec with Matchers with BeforeAndAfterAll with TableDrivenPropertyChecks {

  behavior of "WDL Wom executable"

  val wdlWorkflow: String =
    """
      |task t {
      |  String t1
      |  Int? t2
      |  command { ... }
      |}
      |
      |workflow w {
      |  File w1
      |  String? w2
      |  
      |  scatter(i in range(5)) {
      |    call t as u
      |  }
      |  
      |  call t
      |}
    """.stripMargin

  val namespace = WdlNamespace.loadUsingSource(wdlWorkflow, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
  val graph = namespace.workflow.toWomWorkflowDefinition(isASubworkflow = false)
    .valueOr(errors => fail(s"Failed to build a wom definition: ${errors.toList.mkString(", ")}"))
    .graph

  val w1OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.w1").getOrElse(fail("Failed to find an input node for w1")).singleOutputPort
  val w2OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.w2").getOrElse(fail("Failed to find an input node for w2")).singleOutputPort
  val t1OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.t.t1").getOrElse(fail("Failed to find an input node for t1")).singleOutputPort
  val t2OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.t.t2").getOrElse(fail("Failed to find an input node for t2")).singleOutputPort
  val u1OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.u.t1").getOrElse(fail("Failed to find an input node for u1")).singleOutputPort
  val u2OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.u.t2").getOrElse(fail("Failed to find an input node for u2")).singleOutputPort

  def validate(inputFile: String): Checked[ResolvedExecutableInputs] = {
    namespace.toWomExecutable(Option(inputFile), NoIoFunctionSet, strictValidation = true) match {
      case Left(errors) => Left(errors)
      case Right(e) => e.resolvedExecutableInputs.validNelCheck
    }
  }

  val validations: TableFor3[String, String, Checked[ResolvedExecutableInputs]] = Table(
    ("test name", "inputs JSON", "input set"),
    (
      "Workflow and task inputs 1",
      """
        |{
        |  "w.w1": "my_file.txt",
        |  "w.t.t1": "helloT",
        |  "w.u.t1": "helloU"
        |}
      """.stripMargin,
      Map[OutputPort, ResolvedExecutableInput](
        w1OutputPort -> Coproduct[ResolvedExecutableInput](WomSingleFile("my_file.txt"): WomValue),
        w2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue.none(WomStringType): WomValue),
        t1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloT"): WomValue),
        t2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue.none(WomIntegerType): WomValue),
        u1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloU"): WomValue),
        u2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue.none(WomIntegerType): WomValue)
      ).validNelCheck
    ),
    (
      "Workflow and task inputs 2",
      """
        |{
        |  "w.w1": "my_file.txt",
        |  "w.w2": "inputString",
        |  "w.t.t1": "helloT",
        |  "w.t.t2": 5,
        |  "w.u.t1": "helloU",
        |  "w.u.t2": 6
        |}
      """.stripMargin,
      Map[OutputPort, ResolvedExecutableInput](
        w1OutputPort -> Coproduct[ResolvedExecutableInput](WomSingleFile("my_file.txt"): WomValue),
        w2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue(WomString("inputString")): WomValue),
        t1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloT"): WomValue),
        t2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue(WomInteger(5)): WomValue),
        u1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloU"): WomValue),
        u2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue(WomInteger(6)): WomValue)
      ).validNelCheck
    ),
    (
      "Missing workflow and task inputs",
      """
        |{
        |}
      """.stripMargin,
      NonEmptyList.fromListUnsafe(List(
        "Required workflow input 'w.t.t1' not specified",
        "Required workflow input 'w.u.t1' not specified",
        "Required workflow input 'w.w1' not specified"
      )).asLeft[ResolvedExecutableInputs]
    )
  )

  /*
   * Note that we create the graph twice:
    * once in namespace.workflow.toWomWorkflowDefinition for the expectations
    * once in namespace.toWomExecutable to validate the actual inputs
   * So the output ports in the map won't be object matches. We just check the name equality.
   */
  private def validateInexactPortEquality(actual: ResolvedExecutableInputs, expected: ResolvedExecutableInputs) = {
    def mapToPortName(resolvedExecutableInputs: ResolvedExecutableInputs) = resolvedExecutableInputs map {
      case (k, v) => k.identifier.fullyQualifiedName -> v
    }

    mapToPortName(actual) shouldBe mapToPortName(expected)
  }

  forAll(validations) { (name, inputSource, expectation) =>

    it should s"validate $name" in {
      (validate(inputSource), expectation) match {
        case (Left(actualError), Left(expectedError)) => actualError.toList.toSet -- expectedError.toList.toSet should be(Set.empty)
        case (Left(actualError), _) => fail(s"Expected success but got errors: ${actualError.toList.mkString(", ")}")
        case (Right(actualInputs), Right(expectedInputs)) => validateInexactPortEquality(actualInputs, expectedInputs)
        case (_, Left(expectedError)) => fail(s"Expected errors: '${expectedError.toList.mkString(", ")}' but got success ")
      }
    }
  }
}
