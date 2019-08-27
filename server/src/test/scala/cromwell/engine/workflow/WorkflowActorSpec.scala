package cromwell.engine.workflow

import java.time.OffsetDateTime
import java.util.concurrent.atomic.AtomicInteger

import akka.actor.{Actor, ActorRef, ActorSystem, Props}
import akka.testkit.{EventFilter, TestActorRef, TestFSMRef, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell._
import cromwell.backend.{AllBackendInitializationData, JobExecutionMap}
import cromwell.core._
import cromwell.core.path.{DefaultPathBuilder, PathBuilder, PathBuilderFactory}
import cromwell.engine.{EngineFilesystems, EngineWorkflowDescriptor}
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{WorkflowExecutionAbortedResponse, WorkflowExecutionFailedResponse, WorkflowExecutionSucceededResponse}
import cromwell.engine.workflow.lifecycle.finalization.CopyWorkflowLogsActor
import cromwell.engine.workflow.lifecycle.finalization.WorkflowFinalizationActor.{StartFinalizationCommand, WorkflowFinalizationSucceededResponse}
import cromwell.engine.workflow.lifecycle.initialization.WorkflowInitializationActor.{WorkflowInitializationAbortedResponse, WorkflowInitializationFailedResponse}
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor.MaterializeWorkflowDescriptorFailureResponse
import cromwell.engine.workflow.workflowstore.{StartableState, Submitted, WorkflowHeartbeatConfig, WorkflowToStart}
import cromwell.util.SampleWdl.ThreeStep
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._

class WorkflowActorSpec extends CromwellTestKitWordSpec with WorkflowDescriptorBuilderForSpecs with BeforeAndAfter with Eventually {

  // https://doc.akka.io/docs/akka/current/testing.html#expecting-log-messages
  override implicit val actorSystem = ActorSystem(
    "testsystem",
    ConfigFactory.parseString("""akka.loggers = ["akka.testkit.TestEventListener"]""")
  )

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
  val initialJobCtByRootWf = new AtomicInteger()

  before {
    currentWorkflowId = WorkflowId.randomId()

    copyWorkflowLogsProbe = TestProbe()
  }

  private val workflowHeartbeatConfig = WorkflowHeartbeatConfig(ConfigFactory.load())

  private def createWorkflowActor(state: WorkflowActorState, extraPathBuilderFactory: Option[PathBuilderFactory] = None) = {
    val actor = TestFSMRef(
      factory = new WorkflowActorWithTestAddons(
        finalizationProbe = finalizationProbe,
        workflowId = currentWorkflowId,
        startState = Submitted,
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
        jobTokenDispenserActor = TestProbe().ref,
        workflowStoreActor = system.actorOf(Props.empty),
        workflowHeartbeatConfig = workflowHeartbeatConfig,
        totalJobsByRootWf = initialJobCtByRootWf,
        extraPathBuilderFactory = extraPathBuilderFactory
      ),
      supervisor = supervisorProbe.ref)
    actor.setState(stateName = state, stateData = WorkflowActorData(Option(currentLifecycleActor.ref), Option(descriptor),
      AllBackendInitializationData.empty, StateCheckpoint(InitializingWorkflowState), Submitted))
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
      supervisorProbe.expectNoMessage(AwaitAlmostNothing)
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
      supervisorProbe.expectNoMessage(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Execution succeeds" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch actor
      actor ! WorkflowExecutionSucceededResponse(Map.empty, CallOutputs.empty)
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectNoMessage(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "not run Finalization actor if aborted when in WorkflowUnstartedState" in {
      val actor = createWorkflowActor(WorkflowUnstartedState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      finalizationProbe.expectNoMessage(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "not run Finalization actor if aborted when in MaterializingWorkflowDescriptorState" in {
      val actor = createWorkflowActor(MaterializingWorkflowDescriptorState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      finalizationProbe.expectNoMessage(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "copy workflow logs in the event of MaterializeWorkflowDescriptorFailureResponse" in {
      val actor = createWorkflowActor(MaterializingWorkflowDescriptorState)
      deathwatch watch actor

      copyWorkflowLogsProbe.expectNoMessage(AwaitAlmostNothing)
      actor ! MaterializeWorkflowDescriptorFailureResponse(new Exception("Intentionally failing workflow materialization to test log copying"))
      copyWorkflowLogsProbe.expectMsg(CopyWorkflowLogsActor.Copy(currentWorkflowId, mockDir))

      finalizationProbe.expectNoMessage(AwaitAlmostNothing)
      deathwatch.expectTerminated(actor)
    }

    "log an error when a path builder factory initialization fails" in {
      EventFilter.error(start = "Failed to copy workflow log", pattern = ".*Failing as requested.*", occurrences = 1).intercept {
        val _ = createWorkflowActor(WorkflowSucceededState, Option(new FailingPathBuilderFactory()))
      }
    }

    "log an error when a path builder factory initialization throws" in {
      EventFilter.error(start = "Failed to copy workflow log", pattern = ".*Throwing as requested.*", occurrences = 1).intercept {
        val _ = createWorkflowActor(WorkflowSucceededState, Option(new ThrowingPathBuilderFactory()))
      }
    }
  }
}

class FailingPathBuilderFactory() extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    Future(throw new Exception("Failing as requested"))
  }
}

class ThrowingPathBuilderFactory() extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    throw new Exception("Throwing as requested")
  }
}

class WorkflowActorWithTestAddons(val finalizationProbe: TestProbe,
                                  workflowId: WorkflowId,
                                  startState: StartableState,
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
                                  jobTokenDispenserActor: ActorRef,
                                  workflowStoreActor: ActorRef,
                                  workflowHeartbeatConfig: WorkflowHeartbeatConfig,
                                  totalJobsByRootWf: AtomicInteger,
                                  extraPathBuilderFactory: Option[PathBuilderFactory]) extends WorkflowActor(
  workflowToStart = WorkflowToStart(id = workflowId,
    submissionTime = OffsetDateTime.now,
    state = startState,
    sources = workflowSources,
    hogGroup = HogGroup("foo")),
  conf = conf,
  ioActor = ioActor,
  serviceRegistryActor = serviceRegistryActor,
  workflowLogCopyRouter = workflowLogCopyRouter,
  jobStoreActor = jobStoreActor,
  subWorkflowStoreActor = subWorkflowStoreActor,
  callCacheReadActor = callCacheReadActor,
  callCacheWriteActor = callCacheWriteActor,
  dockerHashActor = dockerHashActor,
  jobTokenDispenserActor = jobTokenDispenserActor,
  backendSingletonCollection = BackendSingletonCollection(Map.empty),
  workflowStoreActor = workflowStoreActor,
  serverMode = true,
  workflowHeartbeatConfig = workflowHeartbeatConfig,
  totalJobsByRootWf = totalJobsByRootWf,
  fileHashCacheActor = None,
  blacklistCache = None) {

  override val pathBuilderFactories: List[PathBuilderFactory] = extraPathBuilderFactory match {
    case Some(pbf) => EngineFilesystems.configuredPathBuilderFactories :+ pbf
    case None => EngineFilesystems.configuredPathBuilderFactories
  }

  override def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor, jobExecutionMap: JobExecutionMap, worfklowOutputs: CallOutputs) = finalizationProbe.ref
}
