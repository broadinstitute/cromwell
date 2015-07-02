package cromwell

import java.util.UUID

import akka.pattern.ask
import akka.testkit._
import cromwell.binding._
import cromwell.binding.values.WdlString
import cromwell.engine._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DummyDataAccess
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld.Addressee

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps


class SimpleWorkflowActorSpec extends CromwellTestkitSpec("SimpleWorkflowActorSpec") {

  private def buildWorkflowFSMRef(sampleWdl: SampleWdl, rawInputsOverride: Option[WorkflowRawInputs] = None):
  TestFSMRef[WorkflowState, WorkflowFailure, WorkflowActor] = {

    val namespace = WdlNamespace.load(sampleWdl.wdlSource)
    val rawInputs = rawInputsOverride.getOrElse(sampleWdl.rawInputs)
    val coercedInputs = namespace.coerceRawInputs(rawInputs).get
    val descriptor = WorkflowDescriptor(UUID.randomUUID(), namespace, sampleWdl.wdlSource, sampleWdl.wdlJson, coercedInputs)
    TestFSMRef(new WorkflowActor(descriptor, new LocalBackend, DummyDataAccess()))
  }

  override def afterAll() {
    shutdown()
  }

  val TestExecutionTimeout = 5000 milliseconds

  "A WorkflowActor" should {
    "start" in {
      val fsm = buildWorkflowFSMRef(SampleWdl.HelloWorld)
      assert(fsm.stateName == WorkflowSubmitted)
      startingCallsFilter("hello.hello") {
        fsm ! Start
        within(TestExecutionTimeout) {
          awaitCond(fsm.stateName == WorkflowRunning)
          awaitCond(fsm.stateName == WorkflowSucceeded)
          val outputName = "hello.hello.salutation"
          val outputs = Await.result(fsm.ask(GetOutputs).mapTo[WorkflowOutputs], 5 seconds)
          val salutation = outputs.getOrElse(outputName, throw new RuntimeException(s"Output '$outputName' not found."))
          val actualOutput = salutation.asInstanceOf[WdlString].value.trim
          actualOutput shouldEqual "Hello world!"
        }
      }
    }

    "fail to construct with missing inputs" in {
      intercept[UnsatisfiedInputsException] {
        buildWorkflowFSMRef(SampleWdl.HelloWorld, rawInputsOverride = Some(Map.empty))
      }
    }

    "fail to construct with inputs of the wrong type" in {
      intercept[UnsatisfiedInputsException] {
        buildWorkflowFSMRef(SampleWdl.HelloWorld, rawInputsOverride = Some(Map(Addressee -> 3)))
      }
    }

    "fail when a call fails" in {
      val fsm = buildWorkflowFSMRef(SampleWdl.GoodbyeWorld)
      assert(fsm.stateName == WorkflowSubmitted)
      startingCallsFilter("goodbye.goodbye") {
        waitForPattern("persisting status of calls goodbye.goodbye to Starting.") {
          waitForPattern("persisting status of calls goodbye.goodbye to Running.") {
            waitForPattern("persisting status of calls goodbye.goodbye to Failed.") {
              waitForPattern("WorkflowActor .+ transitioning from Running to Failed\\.") {
                fsm ! Start
                awaitCond(fsm.stateName == WorkflowRunning)
                awaitCond(fsm.stateName == WorkflowFailed)
              }
            }
          }
        }
      }
    }
  }
}
