package cromwell

import akka.testkit._
import cromwell.binding.values.WdlString
import cromwell.util.SampleWdl

import scala.language.postfixOps

class OptionalParamWorkflowSpec extends CromwellTestkitSpec("OptionalParamWorkflowSpec") {
  "A workflow with optional parameters" should {
    val event = EventFilter.info(pattern = s"starting calls: optional.hello, optional.hello2, optional.hello_person", occurrences = 1)
    "accept optional parameters" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OptionalParamWorkflow,
        event,
        expectedOutputs = Map(
          "optional.hello.greeting" -> WdlString("hello john\n"),
          "optional.hello2.greeting" -> WdlString("hello \n"),
          "optional.hello_person.greeting" -> WdlString("hello world\n")
        )
      )
    }
  }
}
