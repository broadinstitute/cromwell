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
        | # enable this when InputRecordSchemas are enabled      
        | #w9:
        | #  type:
        | #    name: w9
        | #    type: record
        | #    fields:
        | #    - name: w9a
        | #      type: record
        | #      fields: 
        | #      - name: w9aa
        | #        type: string  
        | w10: Directory
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

  lazy val graph = cwlWorkflow.womDefinition(AcceptAllRequirements, Vector.empty) match {
    case Left(errors) => fail(s"Failed to build wom definition: ${errors.toList.mkString(", ")}")
    case Right(womDef) => womDef.graph
  }

  def getOutputPort(n: Int) = graph.inputNodes.find(_.localName == s"w$n").getOrElse(fail(s"Failed to find an input node for w$n")).singleOutputPort
  
  lazy val w0OutputPort = getOutputPort(0)
  lazy val w1OutputPort = getOutputPort(1)
  lazy val w2OutputPort = getOutputPort(2)
  lazy val w3OutputPort = getOutputPort(3)
  lazy val w4OutputPort = getOutputPort(4)
  lazy val w5OutputPort = getOutputPort(5)
  lazy val w6OutputPort = getOutputPort(6)
  lazy val w7OutputPort = getOutputPort(7)
  lazy val w8OutputPort = getOutputPort(8)
//  lazy val w9OutputPort = getOutputPort(9)
  lazy val w10OutputPort = getOutputPort(10)
  
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
          secondaryFiles:
            - class: File
              path: secondaryFile.txt
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
        w9: 
          w9a:
            w9aa: "hello"
        w10:
          class: Directory
          location: "directory_location"
          listing:
            - class: Directory
              path: "innerDirectory"
            - class: File
              location: "innerFile"
      """.stripMargin

    val validInputs = validate(inputFile).map {
      case (port, resolvedInput) => (port.name, resolvedInput)
    }

    // w0 has no input value in the input file, so it should fallback to the default value
    // TODO WOM: when we have string value for wom expression, check that it's "hi !"
    validInputs(w0OutputPort.name).select[WomExpression].get.sourceString shouldBe "\"hi w0 !\""
    validInputs(w1OutputPort.name) shouldBe
      Coproduct[ResolvedExecutableInput](WomMaybePopulatedFile(
        valueOption = Option("my_file.txt"),
        secondaryFiles = List(WomMaybePopulatedFile("secondaryFile.txt"))
      ): WomValue)
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
    // Enable this when InputRecordSchema is enabled
//    validInputs(w9OutputPort.name) shouldBe Coproduct[ResolvedExecutableInput](
//      WomObject(
//        Map(
//          "w9a" -> WomObject(Map("w9aa" -> WomString("hello")))
//        )
//      ): WomValue
//    )
    validInputs(w10OutputPort.name) shouldBe
      Coproduct[ResolvedExecutableInput](WomMaybeListedDirectory(
        valueOption = Option("directory_location"),
        listingOption = Option(
          List(
            WomMaybeListedDirectory("innerDirectory"),
            WomMaybePopulatedFile("innerFile")
          )
        )
      ): WomValue)
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
        "Required workflow input 'w7' not specified",
        "Required workflow input 'w8' not specified",
        "Required workflow input 'w10' not specified"
      )
    }
  }
}
