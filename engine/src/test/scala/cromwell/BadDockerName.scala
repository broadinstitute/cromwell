package cromwell

import akka.testkit._
import cromwell.core.Tags.DockerTest
import cromwell.core.WorkflowFailed
import cromwell.util.SampleWdl

import scala.language.postfixOps

class BadDockerName extends CromwellTestkitSpec {
  "A task which has a bad docker image name" should {
    "fail properly" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"transition from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1),
        runtime = """
                    |runtime {
                    |  docker: "/fauxbuntu:nosuchversion"
                    |}
                  """.stripMargin,
        expectedOutputs = Map(),
        terminalState = WorkflowFailed
      )
    }
  }
}
