package cromwell

import akka.testkit.EventFilter
import cromwell.core.Tags.DockerTest
import cromwell.util.SampleWdl
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.values.{WdlArray, WdlFile, WdlString}


class InputLocalizationWorkflowSpec extends CromwellTestkitSpec {
  "a workflow running on a SharedFileSystem" should {

    val expectedOutputs = Map(
      "wf.fromDifferentDirectories.ls" -> WdlString("1"),
      "wf.echo_int.out" -> WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("out"), WdlFile("out"))),
      "wf.fromSameDirectory.ls" -> WdlString("2")
    )

    "ensure task inputs isolation" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.InputIsolationWdl,
        eventFilter = EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1),
        expectedOutputs = expectedOutputs
      )
    }

    "ensure task inputs isolation in a docker container" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.InputIsolationWdl,
        eventFilter = EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1),
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
