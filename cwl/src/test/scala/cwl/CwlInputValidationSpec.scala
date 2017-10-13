package cwl

import better.files.{File => BFile}
import cwl.CwlDecoder.decodeAllCwl
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import shapeless.Coproduct
import wom.expression.WomExpression
import wom.graph.Graph.ResolvedExecutableInput
import wom.graph.GraphNodePort
import wom.values._

class CwlInputValidationSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks with BeforeAndAfterAll {
  behavior of "CWL Wom executable"

  var cwlFile: BFile = _
  
  override def beforeAll(): Unit = {
    cwlFile = BFile.newTemporaryFile().write(
      """
        |cwlVersion: v1.0
        |class: Workflow
        |inputs:
        | w0:
        |   type: string
        |   default: "hi w0 !"
        | w1: File
        | w2:
        |   type: string
        |   default: "hi w2 !"
        | w3: int
        | w4: long
        | w5: double
        | w6: float
        | w7: boolean
        |steps: []    
        |outputs: []
      """.stripMargin
    )
  }
  
  override def afterAll(): Unit = {
    cwlFile.delete()
    ()
  }

  lazy val cwlWorkflow = decodeAllCwl(cwlFile).map {
    _.select[Workflow].get
  }.value.unsafeRunSync.fold(error => throw new RuntimeException(s"broken parse! msg was $error"), identity)

  lazy val graph = cwlWorkflow.womDefinition match {
    case Left(errors) => fail(s"Failed to build wom definition: ${errors.toList.mkString(", ")}")
    case Right(womDef) => womDef.graph.getOrElse(fail("Failed to build wom graph"))
  }

  lazy val w0OutputPort = graph.inputNodes.find(_.localName == "w0").getOrElse(fail("Failed to find an input node for w0")).singleOutputPort
  lazy val w1OutputPort = graph.inputNodes.find(_.localName == "w1").getOrElse(fail("Failed to find an input node for w1")).singleOutputPort
  lazy val w2OutputPort = graph.inputNodes.find(_.localName == "w2").getOrElse(fail("Failed to find an input node for w2")).singleOutputPort
  lazy val w3OutputPort = graph.inputNodes.find(_.localName == "w3").getOrElse(fail("Failed to find an input node for w3")).singleOutputPort
  lazy val w4OutputPort = graph.inputNodes.find(_.localName == "w4").getOrElse(fail("Failed to find an input node for w4")).singleOutputPort
  lazy val w5OutputPort = graph.inputNodes.find(_.localName == "w5").getOrElse(fail("Failed to find an input node for w5")).singleOutputPort
  lazy val w6OutputPort = graph.inputNodes.find(_.localName == "w6").getOrElse(fail("Failed to find an input node for w6")).singleOutputPort
  lazy val w7OutputPort = graph.inputNodes.find(_.localName == "w7").getOrElse(fail("Failed to find an input node for w7")).singleOutputPort
  
  def validate(inputFile: String): Map[GraphNodePort.OutputPort, ResolvedExecutableInput] = {
    cwlWorkflow.womExecutable(Option(inputFile)) match {
      case Left(errors) => fail(s"Failed to build a wom executable: ${errors.toList.mkString(", ")}")
      case Right(executable) => executable.resolvedExecutableInputs
    }
  }

  it should "parse and validate a valid input file" in {
    val inputFile =
      """
        w1:
          class: File
          path: my_file.txt
        w2: hello !
        w3: 3
        w4: 4
        w5: 5.1
        w6: 6.1
        w7: true
      """.stripMargin

    val validInputs = validate(inputFile).map {
      case (port, resolvedInput) => (port.name, resolvedInput)
    }

    // w0 has no input value in the input file, so it should fallback to the default value
    // TODO WOM: when we have string value for wom expression, check that it's "hi !"
    validInputs(w0OutputPort.name).select[WomExpression].isDefined shouldBe true
    validInputs(w1OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WdlFile("my_file.txt"): WdlValue)
    validInputs(w2OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WdlString("hello !"): WdlValue)
    validInputs(w3OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WdlInteger(3): WdlValue)
    validInputs(w4OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WdlInteger(4): WdlValue)
    validInputs(w5OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WdlFloat(5.1F): WdlValue)
    validInputs(w6OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WdlFloat(6.1F): WdlValue)
    validInputs(w7OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WdlBoolean(true): WdlValue)
  }

  it should "not validate when required inputs are missing" in {
    val inputFile =
      """
        w2: hello !
      """.stripMargin

    cwlWorkflow.womExecutable(Option(inputFile)) match {
      case Right(booh) => fail(s"Expected failed validation but got valid input map: $booh")
      case Left(errors) => errors.toList.toSet shouldBe Set(
        "Required workflow input 'w1' not specified",
        "Required workflow input 'w3' not specified",
        "Required workflow input 'w4' not specified",
        "Required workflow input 'w5' not specified",
        "Required workflow input 'w6' not specified",
        "Required workflow input 'w7' not specified"
      )
    }
  }
}
