package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import wdl4s.values.WdlInteger
import cromwell.util.SampleWdl

import scala.language.postfixOps

object ThreeStepActorSpec {
  val CannedExpectations = Map(
    "three_step.cgrep.count" -> WdlInteger(3),
    "three_step.wc.count" -> WdlInteger(6)
  )
  val EventMessage = s"starting calls: three_step.cgrep, three_step.wc"
}

class ThreeStepActorSpec extends CromwellTestkitSpec {
  import ThreeStepActorSpec._

  "A three step workflow" should {
    "best get to (three) steppin'" ignore {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CannedThreeStep,
        EventFilter.info(pattern = EventMessage, occurrences = 1),
        expectedOutputs = CannedExpectations)
    }
  }

  "A Dockerized three step workflow" should {
    "best get to Dockerized (three) steppin'" taggedAs DockerTest ignore {
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

    "pass canned files properly" taggedAs DockerTest ignore {
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
