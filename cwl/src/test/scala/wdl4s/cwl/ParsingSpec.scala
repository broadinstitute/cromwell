package cwl

import org.scalacheck.Properties
import io.circe.parser._
import io.circe.refined._
import io.circe.literal._
import CwlCodecs._

class WorkflowParsingSpec extends Properties("Workflow Json Parser") {

  def workflowJson(classValue: String) = s"""{"class":"$classValue", "inputs":[], "outputs":[], "steps":[]}"""

  property("accepts the Workflow argument") =
    decode[Workflow](workflowJson("Workflow")).isRight

  property("doesn't parse w/something other than Workflow as class") =
    decode[Workflow](workflowJson("wrong")).isLeft

  property("doesn't parse when class argument is missing") =
    decode[Workflow]("""{"inputs":[], "outputs":[], "steps":[]}""").isLeft
}

class CommandLineToolParsingSpec extends Properties("CommandLineTool Json Parser") {

  def commandLineToolJson(classValue: String) = s"""{"class":"$classValue", "inputs":[], "outputs":[]}"""

  property("accepts the CommandLineTool argument for class") =
    decode[CommandLineTool](commandLineToolJson("CommandLineTool")).isRight

  property("doesn't parse w/something other than") =
    decode[CommandLineTool](commandLineToolJson("wrong")).isLeft

  property("doesn't parse when class argument is missing") =
    decode[CommandLineTool]("""{"inputs":[], "outputs":[]}""").isLeft
}
