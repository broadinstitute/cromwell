package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import wdl4s.types.{WdlArrayType, WdlStringType, WdlFileType}
import wdl4s.{WorkflowInput, NamespaceWithWorkflow}
import wdl4s.values.{WdlArray, WdlInteger, WdlString}
import cromwell.engine.backend.BackendType
import cromwell.util.SampleWdl

import scala.language.postfixOps

class WriteTsvSpec extends CromwellTestkitSpec {
  val row1 = WdlArray(WdlArrayType(WdlStringType), Seq("a", "b").map(WdlString))
  val row2 = WdlArray(WdlArrayType(WdlStringType), Seq("c", "d").map(WdlString))
  val tsv = WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), Seq(row1, row2))
  val outputs = Map(
    "write_lines.a2f.out0" -> tsv,
    "write_lines.a2f.out1" -> tsv
  )

  "A task that calls write_tsv() called in the command section or declaration" should {
    "run properly" ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WriteTsvWorkflow,
        eventFilter = EventFilter.info(pattern = s"starting calls: write_lines.a2f", occurrences = 1),
        expectedOutputs = outputs
      )
    }
    "run properly in a Docker environment" taggedAs DockerTest ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WriteTsvWorkflow,
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
