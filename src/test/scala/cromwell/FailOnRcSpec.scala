package cromwell

import akka.testkit.EventFilter
import cromwell.engine.WorkflowFailed
import cromwell.util.SampleWdl

class FailOnRcSpec extends CromwellTestkitSpec("FailOnRcSpec") {
  "A workflow with tasks that produce non-zero result codes" should {
    "have correct contents in stdout/stderr files for a call that implicitly fails on rc" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.FailOnRc,
        eventFilter = EventFilter.info(pattern = s"persisting status of w.A to Failed", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq("")),
        terminalState = WorkflowFailed
      )
    }

    "have correct contents in stdout/stderr files for a call that explicitly mentions fail on rc" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.FailOnRc,
        runtime = "runtime {failOnRc: true}",
        eventFilter = EventFilter.info(pattern = s"persisting status of w.A to Failed", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq("")),
        terminalState = WorkflowFailed
      )
    }

    "have correct contents in stdout/stderr files for a call that does not fail on rc" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.FailOnRc,
        runtime = "runtime {failOnRc: false}",
        eventFilter = EventFilter.info(pattern = s"persisting status of w.B to Done", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n"), "w.B" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq(""), "w.B" -> Seq(""))
      )
    }
  }
}
