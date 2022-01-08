package cromwell

import cromwell.core.WorkflowFailed
import cromwell.util.SampleWdl


// TODO: These tests are (and were) somewhat unsatisfactory. They'd be much better if we use TestFSMRefs and TestProbes to simulate job completions against the WorkflowActor and make sure it only completes the workflow at the appropriate time.
class WorkflowFailSlowSpec extends CromwellTestKitWordSpec {
  val FailFastOptions: String =
    """
      |{
      |  "workflow_failure_mode": "NoNewCalls"
      |}
    """.stripMargin

  "A workflow containing a failing task" should {
    "not complete any other tasks and ultimately fail, for NoNewCalls" in {
      val outputs = runWdl(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        workflowOptions = FailFastOptions,
        terminalState = WorkflowFailed,
        testActorName = "TestCromwellRootActor-not-complete",
      )
      outputs.size should be(0)
    }
  }

  "A workflow containing a failing task" should {
    "behave like NoNewCalls, if no workflowFailureMode is set" in {
      val outputs = runWdl(
        sampleWdl = SampleWdl.WorkflowFailSlow,
        terminalState = WorkflowFailed,
        testActorName = "TestCromwellRootActor-workflowFailureMode",
      )
      outputs.size should be(0)
    }
  }
}
