package cromwell.engine.workflow

import java.nio.file.Paths

import akka.actor.{Actor, ActorRef}
import akka.testkit.{TestActorRef, TestFSMRef, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.{AllBackendInitializationData, JobExecutionMap}
import cromwell.core._
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.{CopyWorkflowLogsActor, EngineLifecycleActorAbortCommand}
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.MaterializeWorkflowDescriptorFailureResponse
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{WorkflowInitializationAbortedResponse, WorkflowInitializationFailedResponse}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{WorkflowExecutionAbortedResponse, WorkflowExecutionFailedResponse, WorkflowExecutionSucceededResponse}
import cromwell.util.SampleWdl.ThreeStep
import cromwell._
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._

class WorkflowActorSpec extends CromwellTestKitSpec with WorkflowDescriptorBuilder with BeforeAndAfter with Eventually {
  override implicit val actorSystem = system

  val mockServiceRegistryActor = TestActorRef(new Actor {
    override def receive = {
      case _ => // No action
    }
  })

  val mockDir = Paths.get("/where/to/copy/wf/logs")
  val mockWorkflowOptions = s"""{ "final_workflow_log_dir" : "$mockDir" }"""

  var currentWorkflowId: WorkflowId = _
  val currentLifecycleActor = TestProbe()
  val wdlSources = ThreeStep.asWorkflowSources(workflowOptions = mockWorkflowOptions)
  val descriptor = createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = wdlSources)
  val supervisorProbe = TestProbe()
  val deathwatch = TestProbe()
  val finalizationProbe = TestProbe()
  val copyWorkflowLogsProbe = TestProbe()
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
        workflowLogCopyRouter = copyWorkflowLogsProbe.ref,
        jobStoreActor = system.actorOf(AlwaysHappyJobStoreActor.props),
        subWorkflowStoreActor = system.actorOf(AlwaysHappySubWorkflowStoreActor.props),
        callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props),
        jobTokenDispenserActor = TestProbe().ref
      ),
      supervisor = supervisorProbe.ref)
    actor.setState(stateName = state, stateData = WorkflowActorData(Option(currentLifecycleActor.ref), Option(descriptor),
      AllBackendInitializationData.empty, StateCheckpoint(InitializingWorkflowState)))
    actor
  }

  implicit val TimeoutDuration = CromwellTestKitSpec.TimeoutDuration

  "WorkflowActor" should {

    "run Finalization actor if Initialization fails" in {
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
      actor ! WorkflowExecutionFailedResponse(Map.empty, new Exception("Execution Failed"))
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
      currentLifecycleActor.expectMsgPF(CromwellTestKitSpec.TimeoutDuration) {
        case EngineLifecycleActorAbortCommand =>
          actor ! WorkflowExecutionAbortedResponse(Map.empty)
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
      actor ! WorkflowExecutionSucceededResponse(Map.empty, Map.empty)
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

    "copy workflow logs in the event of MaterializeWorkflowDescriptorFailureResponse" in {
      val actor = createWorkflowActor(MaterializingWorkflowDescriptorState)
      deathwatch watch actor
      actor ! MaterializeWorkflowDescriptorFailureResponse(new Exception("Intentionally failing workflow materialization to test log copying"))
      finalizationProbe.expectNoMsg(AwaitAlmostNothing)
      copyWorkflowLogsProbe.expectMsg(CopyWorkflowLogsActor.Copy(currentWorkflowId, mockDir))
      deathwatch.expectTerminated(actor)
    }
  }
}

class MockWorkflowActor(val finalizationProbe: TestProbe,
                        workflowId: WorkflowId,
                        startMode: StartMode,
                        workflowSources: WorkflowSourceFilesCollection,
                        conf: Config,
                        serviceRegistryActor: ActorRef,
                        workflowLogCopyRouter: ActorRef,
                        jobStoreActor: ActorRef,
                        subWorkflowStoreActor: ActorRef,
                        callCacheReadActor: ActorRef,
                        jobTokenDispenserActor: ActorRef) extends WorkflowActor(workflowId, startMode, workflowSources, conf, serviceRegistryActor, workflowLogCopyRouter, jobStoreActor, subWorkflowStoreActor, callCacheReadActor, jobTokenDispenserActor, BackendSingletonCollection(Map.empty), serverMode = true) {

  override def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor, jobExecutionMap: JobExecutionMap, worfklowOutputs: CallOutputs) = finalizationProbe.ref
}
