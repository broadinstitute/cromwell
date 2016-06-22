package cromwell

import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class StdoutStderrWorkflowSpec extends CromwellTestkitSpec {
  "A workflow with tasks that produce stdout/stderr" should {
    "have correct contents in stdout/stderr files for a call" in {
      runWdlAndAssertLogs(
        sampleWdl = SampleWdl.StdoutStderr,
        runtime = """runtime { failOnStderr: "false" }""",
        eventFilter = workflowSuccessFilter,
        expectedStdoutContent = "somethingOnStdout\n",
        expectedStderrContent = "somethingOnStderr\n"
      )
    }

    "have correct contents in stdout/stderr files in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertLogs(
        sampleWdl = SampleWdl.StdoutStderr,
        runtime = """runtime { failOnStderr: "false" docker: "ubuntu:latest" }""",
        eventFilter = workflowSuccessFilter,
        expectedStdoutContent = "somethingOnStdout\n",
        expectedStderrContent = "somethingOnStderr\n"
      )
    }
  }
}
