package cromwell

import akka.testkit._
import wdl4s.values.{WdlInteger, WdlFile}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class StringInterpolationWorkflowSpec extends CromwellTestkitSpec {
  "A workflow with a task that uses string interpolation" should {
    "interpolate strings correctly and run" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.StringInterpolation,
        EventFilter.info(pattern = s"Starting calls: echo_wf.echo", occurrences = 1),
        expectedOutputs = Map(
          "echo_wf.echo.outfile" -> WdlFile("foobar.txt"),
          "echo_wf.echo.two" -> WdlInteger(2)
        )
      )
    }
  }
}
