package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import cromwell.binding.values.WdlString
import cromwell.util.SampleWdl

import scala.language.postfixOps

class GlobbingWorkflowSpec extends CromwellTestkitSpec("GlobbingWorkflowSpec") {
  def doTheTest(runtime: String = "") = runWdlAndAssertOutputs(
    sampleWdl = SampleWdl.GlobtasticWorkflow,
    eventFilter = EventFilter.info(pattern = s"starting calls: w.B", occurrences = 1),
    runtime = runtime,
    expectedOutputs = Map(
      "w.B.B_out" -> WdlString(s"a${newline}b${newline}c")
    )
  )

  val newline = System.lineSeparator
  "A workflow with globbed outputs" should {
    "run properly" in doTheTest()
    "run properly in a Docker environment" taggedAs DockerTest in doTheTest("""runtime { docker: "ubuntu:latest" }""")
  }
}
