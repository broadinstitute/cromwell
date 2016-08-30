package cromwell

import akka.testkit._
import cromwell.core.WorkflowFailed
import cromwell.util.SampleWdl

class WorkflowFailSlowSpec extends CromwellTestkitSpec {
  val FailFastOptions =
    """
      |{
      |  "workflow_failure_mode": "NoNewCalls"
      |}
    """.stripMargin

  // Equivalent centaur test currently ignored
  "A workflow containing a failing task" should {
    "not complete any other tasks and ultimately fail, for NoNewCalls" in {
      runWdl(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        workflowOptions = FailFastOptions,
        eventFilter = EventFilter.info(pattern = s"Job wf.E:NA:1 succeeded!", occurrences = 0),
        runtime = "",
        terminalState = WorkflowFailed
      )
    }
  }

  "A workflow containing a failing task" should {
    "behave like NoNewCalls, if no workflowFailureMode is set" in {
      runWdl(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        eventFilter = EventFilter.info(pattern = s"Job wf.E:NA:1 succeeded!", occurrences = 0),
        runtime = "",
        terminalState = WorkflowFailed
      )
    }
  }
}
