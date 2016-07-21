package cromwell

import akka.testkit._
import wdl4s.types.{WdlStringType, WdlArrayType}
import wdl4s.values.{WdlArray, WdlString}
import cromwell.core.Tags.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ReadLinesFunctionSpec extends CromwellTestkitSpec {
  val outputArray = WdlArray(WdlArrayType(WdlStringType), Seq(
    WdlString("java"),
    WdlString("scala"),
    WdlString("c"),
    WdlString("c++"),
    WdlString("python"),
    WdlString("bash")
  ))

  "A workflow with a read_lines() call in it" should {
    "convert an output file to an Array[String]" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ReadLinesFunctionWdl,
        eventFilter = EventFilter.info(pattern =
          "Starting calls: read_lines.cat_to_file:NA:1, read_lines.cat_to_stdout:NA:1", occurrences = 1),
        expectedOutputs = Map(
          "read_lines.cat_to_file.lines" -> outputArray,
          "read_lines.cat_to_stdout.lines" -> outputArray
        )
      )
    }
    "convert an output file to an Array[String] in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ReadLinesFunctionWdl,
        eventFilter = EventFilter.info(pattern =
          "Starting calls: read_lines.cat_to_file:NA:1, read_lines.cat_to_stdout:NA:1", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = Map(
          "read_lines.cat_to_file.lines" -> outputArray,
          "read_lines.cat_to_stdout.lines" -> outputArray
        )
      )
    }
  }
}
