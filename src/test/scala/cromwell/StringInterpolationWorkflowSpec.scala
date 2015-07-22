package cromwell

import java.nio.file.Paths

import akka.testkit._
import cromwell.binding.values.{WdlFile, WdlString}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class StringInterpolationWorkflowSpec extends CromwellTestkitSpec("StringInterpolationWorkflowSpec") {
  "A workflow with a task that uses string interpolation" should {
    "interpolate strings correctly and run" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.StringInterpolation,
        EventFilter.info(pattern = s"starting calls: echo_wf.echo", occurrences = 1),
        expectedOutputs = Map(
          "echo_wf.echo.outfile" -> WdlFile("foobar.txt")
        )
      )
    }
  }
}
