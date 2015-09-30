package cromwell

import akka.testkit.EventFilter
import cromwell.binding.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet}
import cromwell.engine.WorkflowFailed
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks._

class ContinueOnReturnCodeSpec extends CromwellTestkitSpec("ContinueOnReturnCodeSpec") {
  "A workflow with tasks that produce non-zero return codes" should {
    "have correct contents in stdout/stderr files for a call that implicitly continues on return code" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        eventFilter = EventFilter.info(pattern = s"persisting status of w.A to Failed", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq("")),
        terminalState = WorkflowFailed
      )
    }

    "have correct contents in stdout/stderr files for a call that explicitly mentions continue on return code" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: false}",
        eventFilter = EventFilter.info(pattern = s"persisting status of w.A to Failed", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq("")),
        terminalState = WorkflowFailed
      )
    }

    "have correct contents in stdout/stderr files for a call that does not continue on return code flag" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: true}",
        eventFilter = EventFilter.info(pattern = s"persisting status of w.B to Done", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n"), "w.B" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq(""), "w.B" -> Seq(""))
      )
    }

    "have correct contents in stdout/stderr files for a call that does not continue on return code value" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: 123}",
        eventFilter = EventFilter.info(pattern = s"persisting status of w.B to Done", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n"), "w.B" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq(""), "w.B" -> Seq(""))
      )
    }

    "have correct contents in stdout/stderr files for a call that does not continue on return code list" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.ContinueOnReturnCode,
        runtime = "runtime {continueOnReturnCode: [123]}",
        eventFilter = EventFilter.info(pattern = s"persisting status of w.B to Done", occurrences = 1),
        stdout = Map("w.A" -> Seq("321\n"), "w.B" -> Seq("321\n")),
        stderr = Map("w.A" -> Seq(""), "w.B" -> Seq(""))
      )
    }
  }

  "Checking for return codes" should {
    "continue on expected return code flags" in {
      val flagTests = Table(
        ("flag", "returnCode", "expectedContinue"),
        (true, 0, true),
        (true, 1, true),
        (false, 0, true),
        (false, 1, false))

      forAll(flagTests) { (flag, returnCode, expectedContinue) =>
        ContinueOnReturnCodeFlag(flag).continueFor(returnCode) should be(expectedContinue)
      }
    }

    "continue on expected return code sets" in {
      val setTests = Table(
        ("set", "returnCode", "expectedContinue"),
        (Set(0), 0, true),
        (Set(0), 1, false),
        (Set(1), 0, false),
        (Set(1), 1, true),
        (Set(0, 1), 0, true),
        (Set(0, 1), 1, true),
        (Set(0, 1), 2, false))

      forAll(setTests) { (set, returnCode, expectedContinue) =>
        ContinueOnReturnCodeSet(set).continueFor(returnCode) should be(expectedContinue)
      }
    }
  }
}
