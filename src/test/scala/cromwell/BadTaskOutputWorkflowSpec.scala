package cromwell

import akka.testkit._
import cromwell.CromwellSpec.{IntegrationTest, DockerTest}
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

    "fail properly in a unknown Docker environment" taggedAs DockerTest in {
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
