package cromwell

import akka.testkit._
import wdl4s.values.{WdlFile, WdlString}
import cromwell.util.SampleWdl

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

class FilePassingWorkflowSpec extends CromwellTestkitSpec {
  "A workflow that passes files between tasks" should {
    "pass files properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingWorkflow,
        EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "file_passing.a.out" -> WdlFile("out"),
          "file_passing.a.out_interpolation" -> WdlFile("out"),
          "file_passing.a.contents" -> WdlString("foo bar baz"),
          "file_passing.b.out" -> WdlFile("out"),
          "file_passing.b.out_interpolation" -> WdlFile("out"),
          "file_passing.b.contents" -> WdlString("foo bar baz")
        ),
        patienceConfig = PatienceConfig(2.minutes.dilated)
      )
    }
  }
}
