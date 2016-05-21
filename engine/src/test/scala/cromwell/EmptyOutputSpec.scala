package cromwell

import akka.testkit.EventFilter
import cromwell.util.SampleWdl
import wdl4s.values.WdlString

class EmptyOutputSpec extends CromwellTestkitSpec {

  "Workflow Actor" should {
    "Run a call with an empty string as output" in {
      runWdlAndAssertOutputs(SampleWdl.EmptyString,
        EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        SampleWdl.EmptyString.outputMap
      )
    }
  }

}
