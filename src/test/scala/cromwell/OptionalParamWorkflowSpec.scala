package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import cromwell.binding.WdlNamespace
import cromwell.binding.values.{WdlFile, WdlString}
import cromwell.parser.BackendType
import cromwell.util.SampleWdl

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

  "A workflow with an optional parameter that has a prefix inside the tag" should {
    "not include that prefix if no value is specified" in {
      val wf = """
         |task find {
         |  command {
         |    find ${File root} ${"-name " pattern?}
         |  }
         |}
         |
         |workflow wf {
         |  call find
         |}
       """.stripMargin
      val ns = WdlNamespace.load(wf, BackendType.LOCAL)
      val findTask = ns.findTask("find") getOrElse {
        fail("Expected to find task 'find'")
      }

      val instantiateWithoutValue = findTask.command.instantiate(Map("root" -> WdlFile("src"))) getOrElse {
        fail("Expected instantiation to work")
      }
      instantiateWithoutValue shouldEqual "find src"

      val instantiateWithValue = findTask.command.instantiate(Map(
        "root" -> WdlFile("src"),
        "pattern" -> WdlString("*.java")
      )).getOrElse {fail("Expected instantiation to work")}
      instantiateWithValue shouldEqual "find src -name *.java"
    }
  }
}
