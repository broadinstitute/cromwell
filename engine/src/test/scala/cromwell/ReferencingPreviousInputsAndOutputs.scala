package cromwell

import akka.testkit._
import cromwell.util.SampleWdl
import wdl4s.values.WdlFloat

import scala.language.postfixOps

class ReferencingPreviousInputsAndOutputs extends CromwellTestkitSpec {
  "A task with outputs which reference other outputs" should {
    "run without let or hindrance" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ReferencingPreviousInputsAndOutputs,
        eventFilter = EventFilter.info(pattern = s"Starting calls: wf.golden_pie", occurrences = 1),
        expectedOutputs = Map(
          "wf.golden_pie.Au" -> WdlFloat(1.6180339887),
          "wf.golden_pie.doubleAu" -> WdlFloat(1.6180339887 * 2),
          "wf.golden_pie.tauValue" -> WdlFloat(3.1415926 * 2),
          "wf.golden_pie.goldenPie" -> WdlFloat(3.1415926 * 1.6180339887)
        )
      )
    }
  }
}
