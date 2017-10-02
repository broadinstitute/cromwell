package cromwell

import akka.testkit._
import cromwell.core.Tags.PostWomTest
import cromwell.util.SampleWdl
import wdl.values.{WdlFile, WdlString}

import scala.concurrent.duration._

class FilePassingWorkflowSpec extends CromwellTestKitWordSpec {
  "A workflow that passes files between tasks" should {
    // TODO WOM: should be fixed once call input evaluation is properly implemented
    "pass files properly" taggedAs PostWomTest ignore {
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
