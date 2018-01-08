package cwl

import better.files.{File => BFile}
import cwl.CwlDecoder.decodeAllCwl
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import shapeless.Coproduct
import wom.expression.WomExpression
import wom.graph.Graph.ResolvedExecutableInput
import wom.graph.GraphNodePort
import wom.types.{WomArrayType, WomStringType}
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
        | w8:
        |   type:
        |     type: array
        |     items:
        |       type: array
        |       items: string
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

  lazy val graph = cwlWorkflow.womDefinition(AcceptAllRequirements) match {
    case Left(errors) => fail(s"Failed to build wom definition: ${errors.toList.mkString(", ")}")
    case Right(womDef) => womDef.graph
  }

  lazy val w0OutputPort = graph.inputNodes.find(_.localName == "w0").getOrElse(fail("Failed to find an input node for w0")).singleOutputPort
  lazy val w1OutputPort = graph.inputNodes.find(_.localName == "w1").getOrElse(fail("Failed to find an input node for w1")).singleOutputPort
  lazy val w2OutputPort = graph.inputNodes.find(_.localName == "w2").getOrElse(fail("Failed to find an input node for w2")).singleOutputPort
  lazy val w3OutputPort = graph.inputNodes.find(_.localName == "w3").getOrElse(fail("Failed to find an input node for w3")).singleOutputPort
  lazy val w4OutputPort = graph.inputNodes.find(_.localName == "w4").getOrElse(fail("Failed to find an input node for w4")).singleOutputPort
  lazy val w5OutputPort = graph.inputNodes.find(_.localName == "w5").getOrElse(fail("Failed to find an input node for w5")).singleOutputPort
  lazy val w6OutputPort = graph.inputNodes.find(_.localName == "w6").getOrElse(fail("Failed to find an input node for w6")).singleOutputPort
  lazy val w7OutputPort = graph.inputNodes.find(_.localName == "w7").getOrElse(fail("Failed to find an input node for w7")).singleOutputPort
  lazy val w8OutputPort = graph.inputNodes.find(_.localName == "w8").getOrElse(fail("Failed to find an input node for w8")).singleOutputPort
  
  def validate(inputFile: String): Map[GraphNodePort.OutputPort, ResolvedExecutableInput] = {
    cwlWorkflow.womExecutable(AcceptAllRequirements, Option(inputFile)) match {
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
          location: my_file.txt
          path: my_file.txt
          path: my_file.txt
        w2: hello !
        w3: 3
        w4: 4
        w5: 5.1
        w6: 6.1
        w7: true
        w8: [
          ["a", "b"],
          ["c", "d"]
        ] 
          
      """.stripMargin

    val validInputs = validate(inputFile).map {
      case (port, resolvedInput) => (port.name, resolvedInput)
    }

    // w0 has no input value in the input file, so it should fallback to the default value
    // TODO WOM: when we have string value for wom expression, check that it's "hi !"
    validInputs(w0OutputPort.name).select[WomExpression].get.sourceString shouldBe "\"hi w0 !\""
    validInputs(w1OutputPort.name) shouldBe
      Coproduct[ResolvedExecutableInput](WomMaybePopulatedFile("my_file.txt"): WomValue)
    validInputs(w2OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WomString("hello !"): WomValue)
    validInputs(w3OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WomInteger(3): WomValue)
    validInputs(w4OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WomInteger(4): WomValue)
    validInputs(w5OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WomFloat(5.1D): WomValue)
    validInputs(w6OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WomFloat(6.1D): WomValue)
    validInputs(w7OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](WomBoolean(true): WomValue)
    validInputs(w8OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](
      WomArray(
        WomArrayType(WomArrayType(WomStringType)),
        List(
          WomArray(WomArrayType(WomStringType), List(WomString("a"), WomString("b"))),
          WomArray(WomArrayType(WomStringType), List(WomString("c"), WomString("d")))
        )
      ): WomValue
    )
  }

  it should "not validate when required inputs are missing" in {
    val inputFile =
      """
        w2: hello !
      """.stripMargin

    cwlWorkflow.womExecutable(AcceptAllRequirements, Option(inputFile)) match {
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
