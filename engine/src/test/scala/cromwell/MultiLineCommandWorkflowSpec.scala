package cromwell

import akka.testkit._
import cromwell.core.Tags.DockerTest
import wdl4s.values.WdlString
import cromwell.util.SampleWdl

import scala.language.postfixOps

class MultiLineCommandWorkflowSpec extends CromwellTestkitSpec {
  "A workflow that calls a task with a multi-line command" should {
    "honor the newlines in the command" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.MultiLineCommandWorkflowWdl,
        eventFilter = EventFilter.info(pattern = "Starting calls: wf.blah", occurrences = 1),
        expectedOutputs = Map("wf.blah.ab" -> WdlString("ab"))
      )
    }
    "honor the newlines in the command in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.MultiLineCommandWorkflowWdl,
        eventFilter = EventFilter.info(pattern = "Starting calls: wf.blah", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = Map("wf.blah.ab" -> WdlString("ab"))
      )
    }
  }
}
