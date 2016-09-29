package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef}
import akka.testkit.{TestActorRef, TestFSMRef, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.AllBackendInitializationData
import cromwell.core.{ExecutionStore, OutputStore, WorkflowId, WorkflowSourceFiles}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{WorkflowInitializationAbortedResponse, WorkflowInitializationFailedResponse}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{WorkflowExecutionAbortedResponse, WorkflowExecutionFailedResponse, WorkflowExecutionSucceededResponse}
import cromwell.util.SampleWdl.ThreeStep
import cromwell.{AlwaysHappyJobStoreActor, CromwellTestkitSpec, EmptyCallCacheReadActor}
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._

class WorkflowActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder with BeforeAndAfter with Eventually {
  override implicit val actorSystem = system

  val mockServiceRegistryActor = TestActorRef(new Actor {
    override def receive = {
      case _ => // No action
    }
  })

  var currentWorkflowId: WorkflowId = _
  val currentLifecycleActor = TestProbe()
  val wdlSources = ThreeStep.asWorkflowSources()
  val descriptor = createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = wdlSources)
  val supervisorProbe = TestProbe()
  val deathwatch = TestProbe()
  val finalizationProbe = TestProbe()

  val AwaitAlmostNothing = 100.milliseconds

  before {
    currentWorkflowId = WorkflowId.randomId()
  }

  private def createWorkflowActor(state: WorkflowActorState) = {
    val actor = TestFSMRef(
      factory = new MockWorkflowActor(
        finalizationProbe = finalizationProbe,
        workflowId = currentWorkflowId,
        startMode = StartNewWorkflow,
        workflowSources = wdlSources,
        conf = ConfigFactory.load,
        serviceRegistryActor = mockServiceRegistryActor,
        workflowLogCopyRouter = TestProbe().ref,
        jobStoreActor = system.actorOf(AlwaysHappyJobStoreActor.props),
        callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props),
        jobTokenDispenserActor = TestProbe().ref
      ),
      supervisor = supervisorProbe.ref)
    actor.setState(stateName = state, stateData = WorkflowActorData(Option(currentLifecycleActor.ref), Option(descriptor),
      AllBackendInitializationData.empty, StateCheckpoint(InitializingWorkflowState)))
    actor
  }

  implicit val TimeoutDuration = CromwellTestkitSpec.TimeoutDuration

  "WorkflowActor" should {

    "run Finalization actor if Initialization fails" in {
      val workflowId = WorkflowId.randomId()
      val actor = createWorkflowActor(InitializingWorkflowState)
      deathwatch watch actor
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Materialization Failed")))
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectMsgPF(TimeoutDuration) { case x: WorkflowFailedResponse => x.workflowId should be(currentWorkflowId) }
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Initialization is aborted" in {
      val actor = createWorkflowActor(InitializingWorkflowState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      eventually { actor.stateName should be(WorkflowAbortingState) }
      currentLifecycleActor.expectMsgPF(TimeoutDuration) {
        case EngineLifecycleActorAbortCommand => actor ! WorkflowInitializationAbortedResponse
      }
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectNoMsg(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization if Execution fails" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch actor
      actor ! WorkflowExecutionFailedResponse(ExecutionStore.empty, OutputStore.empty, Seq(new Exception("Execution Failed")))
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectMsgPF(TimeoutDuration) { case x: WorkflowFailedResponse => x.workflowId should be(currentWorkflowId) }
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Execution is aborted" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      eventually { actor.stateName should be(WorkflowAbortingState) }
      currentLifecycleActor.expectMsgPF(CromwellTestkitSpec.TimeoutDuration) {
        case EngineLifecycleActorAbortCommand =>
          actor ! WorkflowExecutionAbortedResponse(ExecutionStore.empty, OutputStore.empty)
      }
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectNoMsg(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Execution succeeds" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch actor
      actor ! WorkflowExecutionSucceededResponse(ExecutionStore.empty, OutputStore.empty)
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectNoMsg(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "not run Finalization actor if aborted when in WorkflowUnstartedState" in {
      val actor = createWorkflowActor(WorkflowUnstartedState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      finalizationProbe.expectNoMsg(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "not run Finalization actor if aborted when in MaterializingWorkflowDescriptorState" in {
      val actor = createWorkflowActor(MaterializingWorkflowDescriptorState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      finalizationProbe.expectNoMsg(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }
  }
}

class MockWorkflowActor(val finalizationProbe: TestProbe,
                        workflowId: WorkflowId,
                        startMode: StartMode,
                        workflowSources: WorkflowSourceFiles,
                        conf: Config,
                        serviceRegistryActor: ActorRef,
                        workflowLogCopyRouter: ActorRef,
                        jobStoreActor: ActorRef,
                        callCacheReadActor: ActorRef,
                        jobTokenDispenserActor: ActorRef) extends WorkflowActor(workflowId, startMode, workflowSources, conf, serviceRegistryActor, workflowLogCopyRouter, jobStoreActor, callCacheReadActor, jobTokenDispenserActor, BackendSingletonCollection(Map.empty)) {

  override def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor, executionStore: ExecutionStore, outputStore: OutputStore) = finalizationProbe.ref
}
