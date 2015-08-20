package cromwell

import akka.testkit._
import cromwell.binding.values.WdlString
import cromwell.engine.WorkflowFailed
import cromwell.util.SampleWdl

import scala.language.postfixOps

class BadTaskOutputWorkflowSpec extends CromwellTestkitSpec("BadTaskOutputWorkflowSpec") {
  "A task which fails to output a file which we're expecting it to output" should {
    "fail and result in a failed workflow" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.BadTaskOutputWdl,
        EventFilter.info(pattern = s"starting calls: badExample.bad", occurrences = 1),
        expectedOutputs = Map(),
        terminalState = WorkflowFailed
      )
    }
  }
}
