package cwl

import io.circe.parser._
import CwlCodecs._
import org.scalatest.{FlatSpec, Matchers}

class WorkflowParsingSpec extends FlatSpec with Matchers {

  behavior of "Workflow Json Parser"

  def workflowJson(classValue: String) = s"""{"class":"$classValue", "id": "MyCwlWorkflow", "inputs":[], "outputs":[], "steps":[]}"""

  it should "accept the Workflow argument" in {
     decode[Workflow](workflowJson("Workflow")) match {
       case Right(_) => // great!
       case Left(e) => fail(s"Workflow rejected because: ${e.getMessage}")
     }
  }

  it should "not parse w/something other than Workflow as class" in {
    decode[Workflow](workflowJson("wrong")) match {
      case Left(_) => // great!
      case Right(wf) => fail(s"workflow unexpectedly accepted: $wf")
    }
  }

  it should "not parse when class argument is missing" in {
    decode[Workflow]("""{"inputs":[], "outputs":[], "steps":[]}""") match {
      case Left(_) => // great!
      case Right(wf) => fail(s"workflow unexpectedly accepted: $wf")
    }
  }
}

class CommandLineToolParsingSpec extends FlatSpec with Matchers {

  behavior of "CommandLineTool Json Parser"

  def commandLineToolJson(classValue: String) = s"""{"class":"$classValue", "id": "MyCwlTask", "inputs":[], "outputs":[]}"""

  it should "accept the CommandLineTool argument for class" in {
    decode[CommandLineTool](commandLineToolJson("CommandLineTool")) match {
      case Right(_) => // great!
      case Left(e) => fail(s"CommandLineTool rejected because: $e")
    }
  }

  it should "not parse w/something other than CommandLineTool" in {
    decode[CommandLineTool](commandLineToolJson("wrong")) match {
      case Left(_) => // great!
      case Right(clt) => fail(s"CommandLineTool unexpectedly accepted: $clt")
    }
  }

  it should "doesn't parse when class argument is missing" in {
    decode[CommandLineTool]("""{"inputs":[], "outputs":[]}""") match {
      case Left(_) => // great!
      case Right(clt) => fail(s"CommandLineTool unexpectedly accepted: $clt")
    }
  }
}
