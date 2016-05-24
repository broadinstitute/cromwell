package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import wdl4s.types.{WdlArrayType, WdlStringType, WdlFileType}
import wdl4s.{WorkflowInput, NamespaceWithWorkflow}
import wdl4s.values.{WdlArray, WdlInteger, WdlString}
import cromwell.engine.backend.BackendType
import cromwell.util.SampleWdl

import scala.language.postfixOps

class WriteLinesSpec extends CromwellTestkitSpec {
  val outputs = Map(
    "write_lines.a2f.out" -> WdlArray(WdlArrayType(WdlStringType), Seq("a", "b", "c", "d").map(WdlString))
  )

  "A task that calls write_lines() in the command section" should {
    "run properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WriteLinesWorkflow,
        eventFilter = EventFilter.info(pattern = s"starting calls: write_lines.a2f", occurrences = 1),
        expectedOutputs = outputs
      )
    }
    "run properly in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WriteLinesWorkflow,
        eventFilter = EventFilter.info(pattern = s"starting calls: write_lines.a2f", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = outputs
      )
    }
  }
}
