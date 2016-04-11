package cromwell

import akka.testkit._
import cromwell.engine._
import cromwell.engine.backend.WorkflowDescriptorBuilder
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.workflowactor.WorkflowActorMessages.StartNewWorkflow
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.HelloWorld.Addressee

import scala.concurrent.duration._
import scala.language.postfixOps


class SimpleWorkflowActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system

  private def buildWorkflowFSMRef(sampleWdl: SampleWdl, rawInputsOverride: String):
  TestFSMRef[WorkflowState, WorkflowData, WorkflowActor] = {
    val workflowSources = WorkflowSourceFiles(sampleWdl.wdlSource(), rawInputsOverride, "{}")
    val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = workflowSources)
    TestFSMRef(new WorkflowActor(descriptor))
  }

  val TestExecutionTimeout = 10.seconds.dilated

  "A WorkflowActor" should {
    "start, run, succeed and die" in {
      startingCallsFilter("hello.hello") {
        val fsm = buildWorkflowFSMRef(SampleWdl.HelloWorld, SampleWdl.HelloWorld.wdlJson)
        val probe = TestProbe()
        probe watch fsm
        within(TestExecutionTimeout) {
          waitForPattern("transitioning from Submitted to Running") {
            waitForPattern("transitioning from Running to Succeeded") {
              fsm ! StartNewWorkflow()
            }
          }
        }
        probe.expectTerminated(fsm, 10.seconds.dilated)
      }
    }

    "fail to construct with missing inputs" in {
      intercept[IllegalArgumentException] {
        buildWorkflowFSMRef(SampleWdl.HelloWorld, rawInputsOverride = "{}")
      }
    }

    "fail to construct with inputs of the wrong type" in {
      intercept[IllegalArgumentException] {
        buildWorkflowFSMRef(SampleWdl.HelloWorld, rawInputsOverride = s""" { "$Addressee" : 3} """)
      }
    }

    "fail when a call fails" in {
      startingCallsFilter("goodbye.goodbye") {
        waitForPattern("WorkflowActor .+ transitioning from Submitted to Running\\.") {
          waitForPattern("persisting status of goodbye to Starting.") {
            waitForPattern("persisting status of goodbye to Running.") {
              waitForPattern("persisting status of goodbye to Failed.") {
                val fsm = buildWorkflowFSMRef(SampleWdl.GoodbyeWorld, SampleWdl.GoodbyeWorld.wdlJson)
                fsm ! StartNewWorkflow()
              }
            }
          }
        }
      }
    }

    "gracefully handle malformed WDL" in {
      within(TestExecutionTimeout) {
        val fsm = buildWorkflowFSMRef(SampleWdl.CoercionNotDefined, SampleWdl.CoercionNotDefined.wdlJson)
        waitForPattern("transitioning from Running to Failed") {
          fsm ! StartNewWorkflow()
        }
      }
    }
  }
}
