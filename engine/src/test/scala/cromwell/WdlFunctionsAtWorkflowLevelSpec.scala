package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wom.types.{WdlMapType, WdlStringType}
import wom.values.{WdlMap, WdlString}


class WdlFunctionsAtWorkflowLevelSpec extends CromwellTestKitWordSpec {
  val outputMap = WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
    WdlString("k1") -> WdlString("v1"),
    WdlString("k2") -> WdlString("v2"),
    WdlString("k3") -> WdlString("v3")
  ))

  "A workflow with a read_lines() and read_map() at the workflow level" should {
    "execute those functions properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WdlFunctionsAtWorkflowLevel,
        eventFilter = EventFilter.info(pattern = "Starting calls: w.a", occurrences = 1),
        expectedOutputs = Map(
          "w.a.x" -> WdlString("one two three four five"),
          "w.a.y" -> outputMap
        )
      )
    }
  }
}
