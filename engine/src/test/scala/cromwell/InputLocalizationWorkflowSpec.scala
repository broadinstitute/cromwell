package cromwell

import akka.testkit.EventFilter
import cromwell.CromwellSpec.DockerTest
import wdl4s.types.{WdlFileType, WdlStringType, WdlArrayType}
import wdl4s.values.{WdlArray, WdlFile, WdlString}
import cromwell.util.SampleWdl


class InputLocalizationWorkflowSpec extends CromwellTestkitSpec {
  "a workflow running on a SharedFileSystem" should {

    val expectedOutputs = Map(
      "wf.fromDifferentDirectories.ls" -> WdlString("1"),
      "wf.echo_int.out" -> WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("out"), WdlFile("out"))),
      "wf.fromSameDirectory.ls" -> WdlString("2")
    )

    "ensure task inputs isolation" ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.InputIsolationWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: wf.echo_int", occurrences = 1),
        expectedOutputs = expectedOutputs
      )
    }

    "ensure task inputs isolation in a docker container"  taggedAs DockerTest ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.InputIsolationWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: wf.echo_int", occurrences = 1),
        runtime = """
                    |runtime {
                    |  docker: "ubuntu:latest"
                    |}
                  """.stripMargin,
        expectedOutputs = expectedOutputs
      )
    }
  }
}
