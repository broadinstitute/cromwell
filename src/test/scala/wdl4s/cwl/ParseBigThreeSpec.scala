package wdl4s.cwl

import wdl4s.cwl._

import io.circe.syntax._
import io.circe._
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import org.scalatest._
import io.circe.yaml.parser
import io.circe.Json

import io.circe.syntax._
import io.circe._
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import shapeless._, poly._//, ops.union._, union._
import shapeless.ops.coproduct._
import cats._, implicits._//, instances._

import io.circe.yaml.parser
import io.circe._

import io.circe.refined._

object P extends Poly1 {
  implicit val x = at[Map[String, WorkflowStep]] {
    println _
  }

  implicit val y = at[Array[WorkflowStep]] {
    println _
  }
}

class ParseBigThreeSpec extends FlatSpec with Matchers {
  val namespace = "cwl"

  it should "parse 1st tool" in {
  val firstTool = """
cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  message:
    type: string
    inputBinding:
      position: 1
outputs: []
"""

  decode[CommandLineTool](parser.parse(firstTool).right.get.noSpaces)
    .isRight shouldBe true
  }

  it should "parse first workflow" in {
    val firstWorkflow = """
cwlVersion: v1.0
class: Workflow
s: Hi
inputs:
  inp: File
  ex: string

outputs:
  classout:
    type: File
    outputSource: compile/classfile

steps:
  untar:
    run: tar-param.cwl
    in:
      tarfile: inp
      extractfile: ex
    out: [example_out]
  compile:
    run: arguments.cwl
    in:
      src: untar/example_out
    out: [classfile]
"""
    decodeCwl(firstWorkflow)
      .isRight should be (true)
  }


}
