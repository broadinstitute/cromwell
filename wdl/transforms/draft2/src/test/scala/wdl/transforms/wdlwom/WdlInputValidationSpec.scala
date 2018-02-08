package wdl.transforms.wdlwom

import cats.syntax.either._
import common.Checked
import common.validation.Checked._
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import shapeless.Coproduct
import wom.executable.Executable.ResolvedExecutableInputs
import wom.graph.Graph.ResolvedExecutableInput
import wom.transforms.{WomExecutableMaker, WomWorkflowDefinitionMaker}
import wom.types._
import wom.values._

class WdlInputValidationSpec(implicit executableMaker: WomExecutableMaker[WdlNamespaceWithWorkflow],
                             workflowDefinitionMaker: WomWorkflowDefinitionMaker[WdlWorkflow]) extends FlatSpec with Matchers with BeforeAndAfterAll with TableDrivenPropertyChecks {

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
  val graph = namespace.workflow.toWomWorkflowDefinition
    .valueOr(errors => fail(s"Failed to build a wom definition: ${errors.toList.mkString(", ")}"))
    .graph

  val w1OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.w1").getOrElse(fail("Failed to find an input node for w1")).singleOutputPort
  val w2OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.w2").getOrElse(fail("Failed to find an input node for w2")).singleOutputPort
  val t1OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.t.t1").getOrElse(fail("Failed to find an input node for t1")).singleOutputPort
  val t2OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.t.t2").getOrElse(fail("Failed to find an input node for t2")).singleOutputPort
  val u1OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.u.t1").getOrElse(fail("Failed to find an input node for u1")).singleOutputPort
  val u2OutputPort = graph.externalInputNodes.find(_.fullyQualifiedName == "w.u.t2").getOrElse(fail("Failed to find an input node for u2")).singleOutputPort

  def validate(inputFile: String): Checked[ResolvedExecutableInputs] = {
    namespace.toWomExecutable(Option(inputFile)) match {
      case Left(errors) => Left(errors)
      case Right(e) => e.resolvedExecutableInputs.validNelCheck
    }
  }

  it should "validate workflow inputs" in {
    val validations = Table(
      ("inputFile", "expectedResult"),
      (
        """
          |{
          |  "w.w1": "my_file.txt",
          |  "w.t.t1": "helloT",
          |  "w.u.t1": "helloU"
          |}
        """.stripMargin,
        Map (
          w1OutputPort -> Coproduct[ResolvedExecutableInput](WomSingleFile("my_file.txt"): WomValue),
          w2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue.none(WomStringType): WomValue),
          t1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloT"): WomValue),
          t2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue.none(WomIntegerType): WomValue),
          u1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloU"): WomValue),
          u2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue.none(WomIntegerType): WomValue)
        ).validNelCheck
      ),
      (
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
        Map (
          w1OutputPort -> Coproduct[ResolvedExecutableInput](WomSingleFile("my_file.txt"): WomValue),
          w2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue(WomString("inputString")): WomValue),
          t1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloT"): WomValue),
          t2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue(WomInteger(5)): WomValue),
          u1OutputPort -> Coproduct[ResolvedExecutableInput](WomString("helloU"): WomValue),
          u2OutputPort -> Coproduct[ResolvedExecutableInput](WomOptionalValue(WomInteger(6)): WomValue)
        ).validNelCheck
      ),
      (
        """
          |{
          |}
        """.stripMargin,
        Set(
          "Required workflow input 'w.t.t1' not specified",
          "Required workflow input 'w.u.t1' not specified",
          "Required workflow input 'w.w1' not specified"
        ).asLeft[ResolvedExecutableInputs]
      )
    )

    forAll(validations) { (inputSource, expectation) =>
      // The order in the Nel is not important, so make it a Set to check that all the expected failure messages are here, regardless of their order
      validate(inputSource).leftMap(_.toList.toSet) shouldBe expectation
    }
  }

}
