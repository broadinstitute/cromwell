package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef}
import akka.testkit.{TestActorRef, TestFSMRef, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell._
import cromwell.backend.{AllBackendInitializationData, JobExecutionMap}
import cromwell.core._
import cromwell.core.path.DefaultPathBuilder
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.MaterializeWorkflowDescriptorFailureResponse
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{WorkflowInitializationAbortedResponse, WorkflowInitializationFailedResponse}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{WorkflowExecutionAbortedResponse, WorkflowExecutionFailedResponse, WorkflowExecutionSucceededResponse}
import cromwell.engine.workflow.lifecycle.{CopyWorkflowLogsActor, EngineLifecycleActorAbortCommand}
import cromwell.util.SampleWdl.ThreeStep
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._

class WorkflowActorSpec extends CromwellTestKitWordSpec with WorkflowDescriptorBuilder with BeforeAndAfter with Eventually {
  override implicit val actorSystem = system

  val mockServiceRegistryActor = TestActorRef(new Actor {
    override def receive = {
      case _ => // No action
    }
  })

  val mockDir = DefaultPathBuilder.get("/where/to/copy/wf/logs")
  val mockWorkflowOptions = s"""{ "final_workflow_log_dir" : "$mockDir" }"""

  var currentWorkflowId: WorkflowId = _
  val currentLifecycleActor = TestProbe()
  val workflowSources = ThreeStep.asWorkflowSources(workflowOptions = mockWorkflowOptions)
  val descriptor = createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = workflowSources)
  val supervisorProbe = TestProbe()
  val deathwatch = TestProbe()
  val finalizationProbe = TestProbe()
  var copyWorkflowLogsProbe: TestProbe = _
  val AwaitAlmostNothing = 100.milliseconds

  before {
    currentWorkflowId = WorkflowId.randomId()

    copyWorkflowLogsProbe = TestProbe()
  }

  private def createWorkflowActor(state: WorkflowActorState) = {
    val actor = TestFSMRef(
      factory = new MockWorkflowActor(
        finalizationProbe = finalizationProbe,
        workflowId = currentWorkflowId,
        startMode = StartNewWorkflow,
        workflowSources = workflowSources,
        conf = ConfigFactory.load,
        ioActor = system.actorOf(SimpleIoActor.props),
        serviceRegistryActor = mockServiceRegistryActor,
        workflowLogCopyRouter = copyWorkflowLogsProbe.ref,
        jobStoreActor = system.actorOf(AlwaysHappyJobStoreActor.props),
        subWorkflowStoreActor = system.actorOf(AlwaysHappySubWorkflowStoreActor.props),
        callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props),
        callCacheWriteActor = system.actorOf(EmptyCallCacheWriteActor.props),
        dockerHashActor = system.actorOf(EmptyDockerHashActor.props),
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

      copyWorkflowLogsProbe.expectNoMsg(AwaitAlmostNothing)
      actor ! MaterializeWorkflowDescriptorFailureResponse(new Exception("Intentionally failing workflow materialization to test log copying"))
      copyWorkflowLogsProbe.expectMsg(CopyWorkflowLogsActor.Copy(currentWorkflowId, mockDir))

      finalizationProbe.expectNoMsg(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }
  }
}

class MockWorkflowActor(val finalizationProbe: TestProbe,
                        workflowId: WorkflowId,
                        startMode: StartMode,
                        workflowSources: WorkflowSourceFilesCollection,
                        conf: Config,
                        ioActor: ActorRef,
                        serviceRegistryActor: ActorRef,
                        workflowLogCopyRouter: ActorRef,
                        jobStoreActor: ActorRef,
                        subWorkflowStoreActor: ActorRef,
                        callCacheReadActor: ActorRef,
                        callCacheWriteActor: ActorRef,
                        dockerHashActor: ActorRef,
                        jobTokenDispenserActor: ActorRef) extends WorkflowActor(workflowId, startMode, workflowSources, conf, ioActor, serviceRegistryActor, workflowLogCopyRouter, jobStoreActor, subWorkflowStoreActor, callCacheReadActor, callCacheWriteActor, dockerHashActor, jobTokenDispenserActor, BackendSingletonCollection(Map.empty), serverMode = true) {

  override def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor, jobExecutionMap: JobExecutionMap, worfklowOutputs: CallOutputs) = finalizationProbe.ref
}
