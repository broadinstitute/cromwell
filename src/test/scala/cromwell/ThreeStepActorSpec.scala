package cromwell

import cromwell.CromwellSpec.DockerTest
import cromwell.binding.values.WdlInteger
import cromwell.util.SampleWdl
import akka.testkit._

import scala.language.postfixOps

object ThreeStepActorSpec {
  val CannedExpectations = Map(
    "three_step.cgrep.count" -> WdlInteger(3),
    "three_step.wc.count" -> WdlInteger(6)
  )
  val eventMessage = s"starting calls: three_step.cgrep, three_step.wc"
}

class ThreeStepActorSpec extends CromwellTestkitSpec("ThreeStepActorSpec") {
  import ThreeStepActorSpec._
  "A three step workflow" should {
    "best get to (three) steppin'" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CannedThreeStep,
        EventFilter.info(pattern = eventMessage, occurrences = 1),
        expectedOutputs = CannedExpectations)
    }
  }

  "A Dockerized three step workflow" should {
    val event = EventFilter.info(pattern = s"starting calls: three_step.cgrep, three_step.wc", occurrences = 1)
    "best get to Dockerized (three) steppin'" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.CannedThreeStep,
        EventFilter.info(pattern = eventMessage, occurrences = 1),
        runtime = """
                    |runtime {
                    |  docker: "ubuntu:latest"
                    |}
                  """.stripMargin,
        expectedOutputs = CannedExpectations)
    }

    "pass files properly" taggedAs DockerTest in {
      val event = EventFilter.info(pattern = s"starting calls: three_step.cgrep, three_step.wc", occurrences = 1)
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.FilePassingThreeStep,
        EventFilter.info(pattern = eventMessage, occurrences = 1)
      )
    }
  }
}
