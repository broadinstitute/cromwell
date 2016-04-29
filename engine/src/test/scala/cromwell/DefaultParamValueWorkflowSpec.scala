package cromwell

import akka.testkit._
import wdl4s.values.WdlString
import cromwell.util.SampleWdl

import scala.language.postfixOps

class DefaultParamValueWorkflowSpec extends CromwellTestkitSpec {
  "A task with a parameter that has a default value" should {
    "accept a value for that parameter" ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.DefaultParameterValueWithValueSpecified,
        EventFilter.info(pattern = s"starting calls: default.hello", occurrences = 1),
        expectedOutputs = Map("default.hello.greeting" -> WdlString("hello alice"))
      )
    }
    "accept NO value for that parameter" ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.DefaultParameterValueWithNOValueSpecified,
        EventFilter.info(pattern = s"starting calls: default.hello", occurrences = 1),
        expectedOutputs = Map("default.hello.greeting" -> WdlString("hello default value"))
      )
    }
  }
}
