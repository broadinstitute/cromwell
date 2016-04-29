package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import cromwell.util.SampleWdl

import scala.language.postfixOps

class StdoutStderrWorkflowSpec extends CromwellTestkitSpec {
  "A workflow with tasks that produce stdout/stderr" should {
    "have correct contents in stdout/stderr files for a call" ignore {
      runWdlAndAssertStdoutStderr(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"persisting status of hello to Done", occurrences = 1),
        fqn = "hello.hello",
        index = None,
        stdout = Some(Seq("Hello world!\n")),
        stderr = Some(Seq(""))
      )
    }
    "have correct contents in stdout/stderr files for a workflow" ignore {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"persisting status of hello to Done", occurrences = 1),
        stdout = Map("hello.hello" -> Seq("Hello world!\n")),
        stderr = Map("hello.hello" -> Seq(""))
      )
    }
    "have correct contents in stdout/stderr files in a Docker environment" taggedAs DockerTest ignore {
      runWdlAndAssertStdoutStderr(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"persisting status of hello to Done", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        fqn = "hello.hello",
        index = None,
        stdout = Some(Seq("Hello world!\n")),
        stderr = Some(Seq(""))
      )
    }
  }
}
