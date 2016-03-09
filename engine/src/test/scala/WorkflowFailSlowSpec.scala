package cromwell

import akka.testkit._
import cromwell.engine.WorkflowFailed
import cromwell.util.SampleWdl

class WorkflowFailSlowSpec extends CromwellTestkitSpec {
  "A workflow containing a failing task" should {
    "complete other tasks but ultimately fail" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        eventFilter = EventFilter.info(pattern = s"persisting status of E to Done.", occurrences = 1),
        runtime = "",
        expectedOutputs = Map(),
        terminalState = WorkflowFailed
      )
    }
  }
}
