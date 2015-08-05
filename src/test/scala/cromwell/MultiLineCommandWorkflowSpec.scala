package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import cromwell.binding.values.WdlString
import cromwell.util.SampleWdl

import scala.language.postfixOps

class MultiLineCommandWorkflowSpec extends CromwellTestkitSpec("MultiLineCommandWorkflowSpec") {
  "A workflow that calls a task with a multi-line command" should {
    "honor the newlines in the command" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.MultiLineCommandWorkflowWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: wf.blah", occurrences = 1),
        expectedOutputs = Map("wf.blah.ab" -> WdlString("ab"))
      )
    }
    "honor the newlines in the command in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.MultiLineCommandWorkflowWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: wf.blah", occurrences = 1),
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
