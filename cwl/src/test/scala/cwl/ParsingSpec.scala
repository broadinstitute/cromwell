package cwl

import CwlDecoder._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class WorkflowParsingSpec extends AnyFlatSpec with Matchers {

  behavior of "Workflow Json Parser"

  def workflowJson(classValue: String) = s"""{"class":"$classValue", "cwlVersion": "v1.0", "id": "MyCwlWorkflow", "inputs":[], "outputs":[], "steps":[]}"""

  it should "accept the Workflow argument" in {
     decodeCwlString(workflowJson("Workflow")).value.unsafeRunSync() match {
       case Right(_) => // great!
       case Left(e) => fail(s"Workflow rejected because: ${e.toList.mkString(", ")}")
     }
  }

  it should "not parse w/something other than Workflow as class" in {
    decodeCwlString(workflowJson("wrong")).value.unsafeRunSync() match {
      case Left(_) => // great!
      case Right(wf) => fail(s"workflow unexpectedly accepted: $wf")
    }
  }

  it should "not parse when class argument is missing" in {
    decodeCwlString("""{"cwlVersion": "v1.0", "inputs":[], "outputs":[], "steps":[]}""").value.unsafeRunSync() match {
      case Left(_) => // great!
      case Right(wf) => fail(s"workflow unexpectedly accepted: $wf")
    }
  }
}

class CommandLineToolParsingSpec extends AnyFlatSpec with Matchers {

  behavior of "CommandLineTool Json Parser"

  def commandLineToolJson(classValue: String) = s"""{"class":"$classValue", "cwlVersion":"v1.0", "id": "MyCwlTask", "inputs":[], "outputs":[]}"""

  it should "accept the CommandLineTool argument for class" in {
    decodeCwlString(commandLineToolJson("CommandLineTool")).value.unsafeRunSync() match {
      case Right(_) => // great!
      case Left(e) => fail(s"CommandLineTool rejected because: $e")
    }
  }

  it should "not parse w/something other than CommandLineTool" in {
    decodeCwlString(commandLineToolJson("wrong")).value.unsafeRunSync() match {
      case Left(_) => // great!
      case Right(clt) => fail(s"CommandLineTool unexpectedly accepted: $clt")
    }
  }

  it should "doesn't parse when class argument is missing" in {
    decodeCwlString("""{"inputs":[], "outputs":[]}""").value.unsafeRunSync() match {
      case Left(_) => // great!
      case Right(clt) => fail(s"CommandLineTool unexpectedly accepted: $clt")
    }
  }
}
