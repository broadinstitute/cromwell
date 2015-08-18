package cromwell

import akka.testkit._
import cromwell.binding.types.{WdlStringType, WdlArrayType}
import cromwell.binding.values.{WdlArray, WdlString}
import cromwell.CromwellSpec.DockerTest
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.util.SampleWdl

import scala.language.postfixOps

class StdoutStderrWorkflowSpec extends CromwellTestkitSpec("StdoutStderrWorkflowSpec") {
  "A workflow with tasks that produce stdout/stderr" should {
    "have correct contents in stdout/stderr files for a call" in {
      runWdlAndAssertStdoutStderr(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"persisting status of hello.hello to Done", occurrences = 1),
        fqn = "hello.hello",
        index = None,
        stdout = Some("Hello world!\n"),
        stderr = Some("")
      )
    }
    "have correct contents in stdout/stderr files for a workflow" in {
      runWdlAndAssertWorkflowStdoutStderr(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"persisting status of hello.hello to Done", occurrences = 1),
        stdout = Map(ExecutionDatabaseKey("hello.hello", None) -> "Hello world!\n"),
        stderr = Map(ExecutionDatabaseKey("hello.hello", None) -> "")
      )
    }
    "have correct contents in stdout/stderr files in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertStdoutStderr(
        sampleWdl = SampleWdl.HelloWorld,
        eventFilter = EventFilter.info(pattern = s"persisting status of hello.hello to Done", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        fqn = "hello.hello",
        index = None,
        stdout = Some("Hello world!\n"),
        stderr = Some("")
      )
    }
  }
}
