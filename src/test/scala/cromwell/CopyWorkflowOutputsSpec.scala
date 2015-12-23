package cromwell

import akka.testkit.EventFilter
import cromwell.util.SampleWdl


class CopyWorkflowOutputsSpec extends CromwellTestkitSpec {
  import CromwellTestkitSpec._

  // TODO not a real test, replace with Thibault's test.
  "CopyWorkflowOutputs" should {
    "not explode" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ThreeStep,
        eventFilter = EventFilter.info(pattern = s"starting calls: three_step.cgrep, three_step.wc", occurrences = 1),
        runtime = "",
        workflowOptions = """ { "workflow_outputs_destination": "/tmp/mlc_outputs" } """,
        expectedOutputs = Map(
          "three_step.ps.procs" -> AnyValueIsFine,
          "three_step.cgrep.count" -> AnyValueIsFine,
          "three_step.wc.count" -> AnyValueIsFine
        ),
        allowOtherOutputs = false
      )
    }
  }

}
