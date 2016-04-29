package cromwell

import akka.testkit.EventFilter
import cromwell.util.SampleWdl
import wdl4s.values.WdlString

class EmptyOutputSpec extends CromwellTestkitSpec {

  "Workflow Actor" should {
    "Run a call with an empty string as output" ignore {
      runWdlAndAssertOutputs(SampleWdl.EmptyString,
      EventFilter.info(pattern = "starting calls: hello.hello", occurrences = 1),
        SampleWdl.EmptyString.outputMap
      )
    }
  }

}
