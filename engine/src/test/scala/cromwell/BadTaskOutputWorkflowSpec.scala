package cromwell

import akka.testkit._
import cromwell.CromwellSpec.{DockerTest, IntegrationTest}
import cromwell.core.WorkflowFailed
import cromwell.util.SampleWdl

import scala.language.postfixOps

class BadTaskOutputWorkflowSpec extends CromwellTestkitSpec {
  "A task which fails to output a file which we're expecting it to output" should {
    "fail and result in a failed workflow" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.BadTaskOutputWdl,
        EventFilter.info(pattern = s"transition from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1),
        expectedOutputs = Map(),
        terminalState = WorkflowFailed
      )
    }

    "fail properly in a unknown Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.BadTaskOutputWdl,
        eventFilter = EventFilter.info(pattern = s"transition from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1),
        runtime = """
                    |runtime {
                    |  docker: "ubuntu:latest"
                    |}
                  """.stripMargin,
        expectedOutputs = Map(),
        terminalState = WorkflowFailed
      )
    }
  }
}
