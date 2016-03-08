package cromwell

import akka.testkit._
import wdl4s.types.{WdlMapType, WdlStringType, WdlArrayType}
import wdl4s.values.{WdlMap, WdlArray, WdlString}
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class WdlFunctionsAtWorkflowLevelSpec extends CromwellTestkitSpec {
  val outputMap = WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
    WdlString("k1") -> WdlString("v1"),
    WdlString("k2") -> WdlString("v2"),
    WdlString("k3") -> WdlString("v3")
  ))

  "A workflow with a read_lines() and read_map() at the workflow level" should {
    "execute those functions properly" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.WdlFunctionsAtWorkflowLevel,
        eventFilter = EventFilter.info(pattern = s"starting calls: w.a", occurrences = 1),
        expectedOutputs = Map(
          "w.a.x" -> WdlString("one two three four five"),
          "w.a.y" -> outputMap
        )
      )
    }
  }
}
