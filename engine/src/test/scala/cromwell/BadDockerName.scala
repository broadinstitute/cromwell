package cromwell

import akka.testkit._
import cromwell.CromwellSpec.{IntegrationTest, DockerTest}
import cromwell.engine.WorkflowFailed
import cromwell.util.SampleWdl

import scala.language.postfixOps

class BadDockerName extends CromwellTestkitSpec {
  "A task which has a bad docker image name" should {
    "fail properly" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"transitioning from Running to Failed.", occurrences = 1),
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
