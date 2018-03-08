package cromwell

import akka.testkit.EventFilter
import cromwell.util.SampleWdl
import wom.types.{WomMapType, WomStringType}
import wom.values.{WomMap, WomString}


class EngineIoFunctionsAtWorkflowLevelSpec extends CromwellTestKitWordSpec {
  val outputMap = WomMap(WomMapType(WomStringType, WomStringType), Map(
    WomString("k1") -> WomString("v1"),
    WomString("k2") -> WomString("v2"),
    WomString("k3") -> WomString("v3")
  ))

  "A workflow with a read_lines() and read_map() at the workflow level" should {
    "execute those functions properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WdlFunctionsAtWorkflowLevel,
        eventFilter = EventFilter.info(pattern = "Starting calls: w.a", occurrences = 1),
        expectedOutputs = Map(
          "w.a.x" -> WomString("one two three four five"),
          "w.a.y" -> outputMap
        )
      )
    }
  }
}
