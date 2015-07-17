package cromwell

import akka.testkit._
import cromwell.binding.types.{WdlStringType, WdlArrayType}
import cromwell.binding.values.{WdlArray, WdlString}
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ReadLinesFunctionSpec extends CromwellTestkitSpec("ReadLinesFunctionSpec") {
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
        eventFilter = EventFilter.info(pattern = s"starting calls: read_lines.cat_to_file, read_lines.cat_to_stdout", occurrences = 1),
        expectedOutputs = Map(
          "read_lines.cat_to_file.lines" -> outputArray,
          "read_lines.cat_to_stdout.lines" -> outputArray
        )
      )
    }
    "convert an output file to an Array[String] in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ReadLinesFunctionWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: read_lines.cat_to_file, read_lines.cat_to_stdout", occurrences = 1),
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
