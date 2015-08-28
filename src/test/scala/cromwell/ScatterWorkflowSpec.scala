package cromwell

import akka.testkit._
import cromwell.binding.types.{WdlIntegerType, WdlStringType, WdlArrayType}
import cromwell.binding.values.{WdlInteger, WdlArray, WdlString}
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ScatterWorkflowSpec extends CromwellTestkitSpec("ScatterWorkflowSpec") {
  "A workflow with a scatter in it" should {
    "run" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.SimpleScatterWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: scatter0.outside_scatter", occurrences = 1),
        expectedOutputs = Map(
          "scatter0.outside_scatter.out" -> WdlInteger(8000),
          "scatter0.inside_scatter.out" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(1, 2, 3, 4, 5).map(WdlInteger(_)))
        )
      )
    }
  }
}
