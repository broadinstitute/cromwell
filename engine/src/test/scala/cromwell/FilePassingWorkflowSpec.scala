package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wom.values.{WomFile, WomString}

import scala.concurrent.duration._

class FilePassingWorkflowSpec extends CromwellTestKitWordSpec {
  "A workflow that passes files between tasks" should {
    // TODO WOM: should be fixed once call input evaluation is properly implemented
    "pass files properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingWorkflow,
        EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "file_passing.a.out" -> WomFile("out"),
          "file_passing.a.out_interpolation" -> WomFile("out"),
          "file_passing.a.contents" -> WomString("foo bar baz"),
          "file_passing.b.out" -> WomFile("out"),
          "file_passing.b.out_interpolation" -> WomFile("out"),
          "file_passing.b.contents" -> WomString("foo bar baz")
        ),
        patienceConfig = PatienceConfig(2.minutes.dilated)
      )
    }
  }
}
