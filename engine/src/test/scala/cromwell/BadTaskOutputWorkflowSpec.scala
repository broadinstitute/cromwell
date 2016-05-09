package cromwell

import akka.testkit._
import cromwell.CromwellSpec.{IntegrationTest, DockerTest}
import cromwell.engine.WorkflowFailed
import cromwell.util.SampleWdl

import scala.language.postfixOps

class BadTaskOutputWorkflowSpec extends CromwellTestkitSpec {
  "A task which fails to output a file which we're expecting it to output" should {
    "fail and result in a failed workflow" ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.BadTaskOutputWdl,
        EventFilter.info(pattern = s"transitioning from Running to Failed.", occurrences = 1),
        expectedOutputs = Map(),
        terminalState = WorkflowFailed
      )
    }

    "fail properly in a unknown Docker environment" taggedAs DockerTest ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.BadTaskOutputWdl,
        eventFilter = EventFilter.info(pattern = s"transitioning from Running to Failed.", occurrences = 1),
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
