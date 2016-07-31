package cromwell

import akka.testkit._
import cromwell.core.Tags.DockerTest
import wdl4s.values.WdlInteger
import cromwell.util.SampleWdl

import scala.language.postfixOps

object ThreeStepActorSpec {
  val CannedExpectations = Map(
    "three_step.cgrep.count" -> WdlInteger(3),
    "three_step.wc.count" -> WdlInteger(5)
  )
  val EventMessage = "Workflow complete"
}

class ThreeStepActorSpec extends CromwellTestkitSpec {
  import ThreeStepActorSpec._

  "A three step workflow" should {
    "best get to (three) steppin'" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CannedThreeStep,
        EventFilter.info(pattern = EventMessage, occurrences = 1),
        expectedOutputs = CannedExpectations)
    }
  }

  "A Dockerized three step workflow" should {
    "best get to Dockerized (three) steppin'" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CannedThreeStep,
        EventFilter.info(pattern = EventMessage, occurrences = 1),
        runtime = """
                    |runtime {
                    |  docker: "ubuntu:latest"
                    |}
                  """.stripMargin,
        expectedOutputs = CannedExpectations)
    }

    "pass canned files properly" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CannedFilePassing,
        eventFilter = EventFilter.info(pattern = EventMessage, occurrences = 1),
        expectedOutputs = Map(
          "three_step.wc.count" -> WdlInteger(3),
          "three_step.cgrep.count" -> WdlInteger(1)
        )
      )
    }
  }
}
