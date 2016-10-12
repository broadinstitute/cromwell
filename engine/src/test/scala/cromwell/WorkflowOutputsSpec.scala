package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import cromwell.CromwellTestKitSpec.AnyValueIsFine


class WorkflowOutputsSpec extends CromwellTestKitSpec {
  "Workflow outputs" should {
    "use all outputs if none are specified" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ThreeStep,
        eventFilter = EventFilter.info(pattern = s"is in a terminal state: WorkflowSucceededState", occurrences = 1),
        runtime = "",
        expectedOutputs = Map(
          "three_step_ps_procs" -> AnyValueIsFine,
          "three_step_cgrep_count" -> AnyValueIsFine,
          "three_step_wc_count" -> AnyValueIsFine
        ),
        allowOtherOutputs = false
      )
    }

    "Respect the workflow output section" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ThreeStepWithOutputsSection,
        eventFilter = EventFilter.info(pattern = s"is in a terminal state: WorkflowSucceededState", occurrences = 1),
        runtime = "",
        expectedOutputs = Map(
          "three_step_cgrep_count" -> AnyValueIsFine,
          "three_step_wc_count" -> AnyValueIsFine
        ),
        allowOtherOutputs = false
      )
    }

    "Not list scatter shards" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SimpleScatterWdl,
        eventFilter = EventFilter.info(pattern = s"is in a terminal state: WorkflowSucceededState", occurrences = 1),
        runtime = "",
        expectedOutputs = Map(
          "scatter0_outside_scatter_out" -> AnyValueIsFine,
          "scatter0_inside_scatter_out" -> AnyValueIsFine
        ),
        allowOtherOutputs = false
      )
    }

    "Not list scatter shards, even for wildcards" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SimpleScatterWdlWithOutputs,
        eventFilter = EventFilter.info(pattern = s"is in a terminal state: WorkflowSucceededState", occurrences = 1),
        runtime = "",
        expectedOutputs = Map(
          "scatter0_inside_scatter_out" -> AnyValueIsFine
        ),
        allowOtherOutputs = false
      )
    }
  }
}
