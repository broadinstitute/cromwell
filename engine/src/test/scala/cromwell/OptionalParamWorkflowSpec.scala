package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl
import wdl4s.WdlNamespace
import wdl4s.expression.NoFunctions
import wdl4s.values.{WdlFile, WdlString}

import scala.language.postfixOps

class OptionalParamWorkflowSpec extends CromwellTestkitSpec {
  val outputs = Map(
    "optional.hello.greeting" -> WdlString("hello john"),
    "optional.hello2.greeting" -> WdlString("hello"),
    "optional.hello_person.greeting" -> WdlString("hello world")
  )

  "A workflow with optional parameters" should {
    "accept optional parameters" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OptionalParamWorkflow,
        eventFilter = EventFilter.info(pattern =
          "Starting calls: optional.hello:NA:1, optional.hello2:NA:1, optional.hello_person:NA:1", occurrences = 1),
        expectedOutputs = outputs
      )
    }
    "accept optional parameters in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.OptionalParamWorkflow,
        eventFilter = EventFilter.info(pattern =
          "Starting calls: optional.hello:NA:1, optional.hello2:NA:1, optional.hello_person:NA:1", occurrences = 1),
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
         |  String? pattern
         |  File root
         |  command {
         |    find ${root} ${"-name " + pattern}
         |  }
         |}
         |
         |workflow wf {
         |  call find
         |}
       """.stripMargin
      val ns = WdlNamespace.load(wf)
      val findTask = ns.findTask("find") getOrElse {
        fail("Expected to find task 'find'")
      }

      val instantiateWithoutValue = findTask.instantiateCommand(Map("root" -> WdlFile("src")), NoFunctions) getOrElse {
        fail("Expected instantiation to work")
      }
      instantiateWithoutValue shouldEqual "find src"

      val instantiateWithValue = findTask.instantiateCommand(Map(
        "root" -> WdlFile("src"),
        "pattern" -> WdlString("*.java")
      ), NoFunctions).getOrElse {fail("Expected instantiation to work")}
      instantiateWithValue shouldEqual "find src -name *.java"
    }
  }
}
