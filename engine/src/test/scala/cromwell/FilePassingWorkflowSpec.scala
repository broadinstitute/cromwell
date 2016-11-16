package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wdl4s.values.{WdlFile, WdlString}

import scala.concurrent.duration._

class FilePassingWorkflowSpec extends CromwellTestKitSpec {
  "A workflow that passes files between tasks" should {
    "pass files properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingWorkflow,
        EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "file_passing_a_out" -> WdlFile("out"),
          "file_passing_a_out_interpolation" -> WdlFile("out"),
          "file_passing_a_contents" -> WdlString("foo bar baz"),
          "file_passing_b_out" -> WdlFile("out"),
          "file_passing_b_out_interpolation" -> WdlFile("out"),
          "file_passing_b_contents" -> WdlString("foo bar baz")
        ),
        patienceConfig = PatienceConfig(2.minutes.dilated)
      )
    }
  }
}
