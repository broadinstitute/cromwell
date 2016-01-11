package cromwell

import akka.testkit._
import wdl4s.values.{WdlString, WdlFile}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class MultipleFilesWithSameNameWorkflowSpec extends CromwellTestkitSpec {
  "A workflow with two file inputs that have the same name" should {
    "not clobber one file with the contents of another" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FileClobber,
        EventFilter.info(pattern = s"starting calls: two.x, two.y", occurrences = 1),
        expectedOutputs = Map(
          "two.x.out" -> WdlString("first file.txt"),
          "two.y.out" -> WdlString("second file.txt")
        )
      )
    }
  }
}
