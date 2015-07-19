package cromwell

import akka.testkit._
import cromwell.binding.values.WdlString
import cromwell.util.SampleWdl
import cromwell.CromwellSpec.DockerTest

import scala.language.postfixOps

class OptionalParamWorkflowSpec extends CromwellTestkitSpec("OptionalParamWorkflowSpec") {
  val outputs = Map(
    "optional.hello.greeting" -> WdlString("hello john"),
    "optional.hello2.greeting" -> WdlString("hello "),
    "optional.hello_person.greeting" -> WdlString("hello world")
  )

  "A workflow with optional parameters" should {
    "accept optional parameters" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OptionalParamWorkflow,
        eventFilter = EventFilter.info(pattern = s"starting calls: optional.hello, optional.hello2, optional.hello_person", occurrences = 1),
        expectedOutputs = outputs
      )
    }
    "accept optional parameters in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OptionalParamWorkflow,
        eventFilter = EventFilter.info(pattern = s"starting calls: optional.hello, optional.hello2, optional.hello_person", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = outputs
      )
    }
  }
}
