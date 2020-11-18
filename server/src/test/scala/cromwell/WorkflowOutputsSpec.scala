package cromwell

import cromwell.util.SampleWdl
import cromwell.CromwellTestKitSpec.AnyValueIsFine

class WorkflowOutputsSpec extends CromwellTestKitWordSpec {
  "Workflow outputs" should {
    "use all outputs if none are specified" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ThreeStep,
        expectedOutputs = Map(
          "three_step.ps.procs" -> AnyValueIsFine,
          "three_step.cgrep.count" -> AnyValueIsFine,
          "three_step.wc.count" -> AnyValueIsFine
        ),
        allowOtherOutputs = false,
        testActorName = "TestCromwellRootActor-use-all",
      )
    }

    "Respect the workflow output section" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ThreeStepWithOutputsSection,
        expectedOutputs = Map(
          "three_step.cgrep.count" -> AnyValueIsFine,
          "three_step.wc.count" -> AnyValueIsFine
        ),
        allowOtherOutputs = false,
        testActorName = "TestCromwellRootActor-output",
      )
    }

    "Not list scatter shards" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SimpleScatterWdl,
        expectedOutputs = Map(
          "scatter0.outside_scatter.out" -> AnyValueIsFine,
          "scatter0.inside_scatter.out" -> AnyValueIsFine
        ),
        allowOtherOutputs = false,
        testActorName = "TestCromwellRootActor-shards",
      )
    }

    "Not list scatter shards, even for wildcards" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SimpleScatterWdlWithOutputs,
        expectedOutputs = Map(
          "scatter0.inside_scatter.out" -> AnyValueIsFine
        ),
        allowOtherOutputs = false,
        testActorName = "TestCromwellRootActor-wildcards",
      )
    }
  }
}
