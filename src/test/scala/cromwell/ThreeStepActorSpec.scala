package cromwell

import java.util.UUID

import akka.pattern.ask
import akka.testkit.TestFSMRef
import cromwell.CromwellSpec.DockerTest
import cromwell.binding._
import cromwell.binding.values.WdlInteger
import cromwell.engine._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DummyDataAccess
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor._
import cromwell.util.SampleWdl

import scala.concurrent.duration._
import scala.language.postfixOps


object ThreeStepActorSpec {
  val TestExecutionTimeout = 5000 milliseconds
  val CannedExpectations = Map(
    "three_step.cgrep.count" -> 3,
    "three_step.wc.count" -> 6)
}

class ThreeStepActorSpec extends CromwellTestkitSpec("ThreeStepActorSpec") {
  import ThreeStepActorSpec._

  private def buildFsmWorkflowActor(sampleWdl: SampleWdl, runtime: String) = {
    val namespace = WdlNamespace.load(sampleWdl.wdlSource(runtime))
    // This is a test and is okay with just throwing if coerceRawInputs returns a Failure.
    val coercedInputs = namespace.coerceRawInputs(sampleWdl.rawInputs).get
    val descriptor = WorkflowDescriptor(UUID.randomUUID(), namespace, sampleWdl.wdlSource(runtime), sampleWdl.wdlJson, coercedInputs)
    val dataAccess = new DummyDataAccess()
    TestFSMRef(new WorkflowActor(descriptor, new LocalBackend, dataAccess))
  }

  private def runAndAssertCorrectness(sampleWdl: SampleWdl, runtime: String = "", expectations: Map[String, Int] = Map.empty): Unit = {
    val fsm = buildFsmWorkflowActor(sampleWdl, runtime)
    assert(fsm.stateName == WorkflowSubmitted)
    startingCallsFilter("three_step.cgrep", "three_step.wc") {
      fsm ! Start
      within(TestExecutionTimeout) {
        awaitCond(fsm.stateName == WorkflowRunning)
        awaitCond(fsm.stateName == WorkflowSucceeded)
        val outputs = fsm.ask(GetOutputs).mapTo[WorkflowOutputs].futureValue

        expectations foreach { case (outputFqn, expectedValue) =>
          val wdlValue = outputs.getOrElse(outputFqn, throw new RuntimeException(s"Output $outputFqn not found"))
          val actualValue = wdlValue.asInstanceOf[WdlInteger].value.toInt
          actualValue shouldEqual expectedValue
        }
      }
    }
  }

  "A three step workflow" should {
    "best get to (three) steppin'" in {
      runAndAssertCorrectness(sampleWdl = SampleWdl.CannedThreeStep, expectations = CannedExpectations)
    }
  }

  "A Dockerized three step workflow" should {
    "best get to Dockerized (three) steppin'" taggedAs DockerTest in {
      runAndAssertCorrectness(
        sampleWdl = SampleWdl.CannedThreeStep,
        runtime = """
                    |runtime {
                    |  docker: "ubuntu:latest"
                    |}
                  """.stripMargin,
        expectations = CannedExpectations)
    }

    "pass files properly" taggedAs DockerTest in {
      runAndAssertCorrectness(sampleWdl = SampleWdl.FilePassingThreeStep)
    }
  }
}
