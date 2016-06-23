package cromwell

import akka.testkit.EventFilter
import cromwell.core.WorkflowFailed
import cromwell.util.SampleWdl

class ContinueOnReturnCodeWorkflowSpec extends CromwellTestkitSpec {
  "A workflow with tasks that produce non-zero return codes" should {
    "Fail if the return code is undefined in the continueOnReturnCode runtime attribute and the return code is non zero" in {
      runWdl(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        eventFilter = EventFilter.info(pattern = "transition from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1),
        terminalState = WorkflowFailed
      )
    }

    "Fail if the return code is false in the continueOnReturnCode runtime attribute and the return code is non zero" in {
      runWdl(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: false}",
        eventFilter = EventFilter.info(pattern = "transition from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1),
        terminalState = WorkflowFailed
      )
    }

    "Succeed if the return code is true in the continueOnReturnCode runtime attribute" in {
      runWdl(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: true}",
        eventFilter = EventFilter.info(pattern = "transition from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1)
      )
    }

    "Succeed if the return code is defined in the continueOnReturnCode runtime attribute" in {
      runWdl(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: 123}",
        eventFilter = EventFilter.info(pattern = "transition from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1)
      )
    }

    "Succeed if the return code is present in the continueOnReturnCode runtime attributes list" in {
      runWdl(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: [123]}",
        eventFilter = EventFilter.info(pattern = "transition from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1)
      )
    }
  }
}
