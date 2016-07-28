package cromwell.engine.workflow

import akka.actor.Actor
import akka.testkit.{EventFilter, TestActorRef, TestFSMRef, TestProbe}
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.backend.AllBackendInitializationData
import cromwell.core.{ExecutionStore, OutputStore, WorkflowId}
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{WorkflowInitializationAbortedResponse, WorkflowInitializationFailedResponse}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{WorkflowExecutionAbortedResponse, WorkflowExecutionFailedResponse, WorkflowExecutionSucceededResponse}
import cromwell.jobstore.{JobStoreActor, WriteCountingJobStoreDatabase}
import cromwell.util.SampleWdl.ThreeStep
import org.scalatest.BeforeAndAfter

class WorkflowActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder with BeforeAndAfter {
  override implicit val actorSystem = system

  val mockServiceRegistryActor = TestActorRef(new Actor {
    override def receive = {
      case _ => // No action
    }
  })

  val currentLifecycleActor = TestProbe()
  val wdlSources = ThreeStep.asWorkflowSources()
  val descriptor = createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = wdlSources)

  private def createWorkflowActor(state: WorkflowActorState) = {
    val actor = TestFSMRef(new WorkflowActor(WorkflowId.randomId(),
      StartNewWorkflow,
      wdlSources,
      ConfigFactory.load,
      mockServiceRegistryActor,
      system.actorOf(JobStoreActor.props(WriteCountingJobStoreDatabase.makeNew)),
      TestProbe().ref))
    actor.setState(stateName = state, stateData = WorkflowActorData(Option(currentLifecycleActor.ref), Option(descriptor),
      AllBackendInitializationData.empty, StateCheckpoint(InitializingWorkflowState)))
    actor
  }

  "WorkflowActor" should {

    "run Finalization actor if Initialization fails" in {
      val actor = createWorkflowActor(InitializingWorkflowState)

      within(CromwellTestkitSpec.TimeoutDuration) {
        EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1).intercept {
          actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Materialization Failed")))
        }
      }

      actor.stop()
    }

    "run Finalization actor if Initialization is aborted" in {
      val actor = createWorkflowActor(InitializingWorkflowState)

      within(CromwellTestkitSpec.TimeoutDuration) {
        EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowAbortedState", occurrences = 1).intercept {
          actor ! AbortWorkflowCommand
          currentLifecycleActor.expectMsgPF(CromwellTestkitSpec.TimeoutDuration) {
            case EngineLifecycleActorAbortCommand => actor ! WorkflowInitializationAbortedResponse
          }
        }
      }

      actor.stop()
    }

    "run Finalization if Execution fails" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)

      within(CromwellTestkitSpec.TimeoutDuration) {
        EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1).intercept {
          actor ! WorkflowExecutionFailedResponse(ExecutionStore.empty, OutputStore.empty,
            Seq(new Exception("Execution Failed")))
        }
      }

      actor.stop()
    }

    "run Finalization actor if Execution is aborted" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)

      within(CromwellTestkitSpec.TimeoutDuration) {
        EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowAbortedState", occurrences = 1).intercept {
          actor ! AbortWorkflowCommand
          currentLifecycleActor.expectMsgPF(CromwellTestkitSpec.TimeoutDuration) {
            case EngineLifecycleActorAbortCommand =>
              actor ! WorkflowExecutionAbortedResponse(ExecutionStore.empty, OutputStore.empty)
          }
        }
      }

      actor.stop()
    }

    "run Finalization actor if Execution succeeds" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)

      within(CromwellTestkitSpec.TimeoutDuration) {
        EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1).intercept {
          actor ! WorkflowExecutionSucceededResponse(ExecutionStore.empty, OutputStore.empty)
        }
      }

      actor.stop()
    }

    "not run Finalization actor if aborted when in WorkflowUnstartedState" in {
      val actor = createWorkflowActor(WorkflowUnstartedState)

      within(CromwellTestkitSpec.TimeoutDuration) {
        EventFilter.info(pattern = "transitioning from WorkflowUnstartedState to WorkflowAbortedState", occurrences = 1).intercept {
          actor ! AbortWorkflowCommand
        }
      }

      actor.stop()
    }

    "not run Finalization actor if aborted when in MaterializingWorkflowDescriptorState" in {
      val actor = createWorkflowActor(MaterializingWorkflowDescriptorState)

      within(CromwellTestkitSpec.TimeoutDuration) {
        EventFilter.info(pattern = "transitioning from MaterializingWorkflowDescriptorState to WorkflowAbortedState", occurrences = 1).intercept {
          actor ! AbortWorkflowCommand
        }
      }

      actor.stop()
    }
  }
}
