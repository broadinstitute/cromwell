package cromwell

import akka.testkit._
import cromwell.engine.WorkflowFailed
import cromwell.util.SampleWdl

class WorkflowFailSlowSpec extends CromwellTestkitSpec {
  val FailSlowOptions =
    """
      |{
      |  "workflowFailureMode": "ContinueWhilePossible"
      |}
    """.stripMargin

  val FailFastOptions =
    """
      |{
      |  "workflowFailureMode": "NoNewCalls"
      |}
    """.stripMargin

  "A workflow containing a failing task" should {
    "complete other tasks but ultimately fail, for ContinueWhilePossible" ignore {
      runWdl(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        workflowOptions = FailSlowOptions,
        eventFilter = EventFilter.info(pattern = s"persisting status of E to Done.", occurrences = 1),
        runtime = "",
        terminalState = WorkflowFailed
      )
    }
  }

  "A workflow containing a failing task" should {
    "not complete any other tasks and ultimately fail, for NoNewCalls" ignore {
      runWdl(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        workflowOptions = FailFastOptions,
        eventFilter = EventFilter.info(pattern = s"persisting status of E to Done.", occurrences = 0),
        runtime = "",
        terminalState = WorkflowFailed
      )
    }
  }

  "A workflow containing a failing task" should {
    "behave like NoNewCalls, if no workflowFailureMode is set" ignore {
      runWdl(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        eventFilter = EventFilter.info(pattern = s"persisting status of E to Done.", occurrences = 0),
        runtime = "",
        terminalState = WorkflowFailed
      )
    }
  }
}
