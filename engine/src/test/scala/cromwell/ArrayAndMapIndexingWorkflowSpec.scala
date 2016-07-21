package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wdl4s.values.{WdlInteger, WdlString}

import scala.language.postfixOps

class ArrayAndMapIndexingWorkflowSpec extends CromwellTestkitSpec {
  "A workflow with array/map indexes in expressions" should {
    "run locally" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ArrayAndMapIndexingWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = Map(
          "test.echo_str.o" -> WdlString("bar"),
          "test.echo_int.o" -> WdlInteger(200)
        )
      )
    }
  }
}
