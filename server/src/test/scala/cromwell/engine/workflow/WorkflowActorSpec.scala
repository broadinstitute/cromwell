package cromwell.engine.workflow

import java.time.OffsetDateTime
import java.util.concurrent.atomic.AtomicInteger
import akka.actor.{Actor, ActorRef, ActorSystem, Kill, Props}
import akka.testkit.{EventFilter, TestActorRef, TestFSMRef, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.StrictLogging
import cromwell._
import cromwell.backend.{AllBackendInitializationData, JobExecutionMap}
import cromwell.core._
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder, PathBuilderFactory}
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowActorWorkComplete
import cromwell.engine.workflow.lifecycle.EngineLifecycleActorAbortCommand
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{
  ExecuteWorkflowCommand,
  WorkflowExecutionAbortedResponse,
  WorkflowExecutionFailedResponse,
  WorkflowExecutionSucceededResponse
}
import cromwell.engine.workflow.lifecycle.finalization.{CopyWorkflowLogsActor, WorkflowCallbackActor}
import cromwell.engine.workflow.lifecycle.finalization.WorkflowFinalizationActor.{
  StartFinalizationCommand,
  WorkflowFinalizationFailedResponse,
  WorkflowFinalizationSucceededResponse
}
import cromwell.engine.workflow.lifecycle.initialization.WorkflowInitializationActor.{
  StartInitializationCommand,
  WorkflowInitializationAbortedResponse,
  WorkflowInitializationFailedResponse,
  WorkflowInitializationSucceededResponse
}
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor.MaterializeWorkflowDescriptorFailureResponse
import cromwell.engine.workflow.workflowstore.{StartableState, Submitted, WorkflowHeartbeatConfig, WorkflowToStart}
import cromwell.engine.{EngineFilesystems, EngineWorkflowDescriptor}
import cromwell.services.metadata.MetadataService.{MetadataWriteSuccess, PutMetadataActionAndRespond}
import cromwell.util.SampleWdl.ThreeStep
import cromwell.util.WomMocks
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually
import wom.values.WomString

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

class WorkflowActorSpec
    extends CromwellTestKitWordSpec
    with WorkflowDescriptorBuilderForSpecs
    with BeforeAndAfter
    with Eventually
    with StrictLogging {

  override protected lazy val actorSystemConfig: Config =
    ConfigFactory.parseString("""akka.loggers = ["akka.testkit.TestEventListener"]""")
  // https://doc.akka.io/docs/akka/current/testing.html#expecting-log-messages

  implicit override lazy val actorSystem: ActorSystem = system

  val mockServiceRegistryActor: TestActorRef[Actor] =
    TestActorRef(
      new Actor {
        override def receive: Receive = {
          case PutMetadataActionAndRespond(events, replyTo, _) =>
            replyTo ! MetadataWriteSuccess(events)
          case _ => // No action
        }
      },
      "mockServiceRegistryActor"
    )

  val mockDir: Path = DefaultPathBuilder.get("/where/to/copy/wf/logs")
  val mockUri: String = "http://example.com"
  val mockWorkflowOptions = s"""{ "final_workflow_log_dir" : "$mockDir", "workflow_callback_uri": "$mockUri" }"""

  var currentWorkflowId: WorkflowId = _
  val currentLifecycleActor: TestProbe = TestProbe("currentLifecycleActor")
  val workflowSources: WorkflowSourceFilesCollection =
    ThreeStep.asWorkflowSources(workflowOptions = mockWorkflowOptions)
  lazy val descriptor: EngineWorkflowDescriptor =
    createMaterializedEngineWorkflowDescriptor(WorkflowId.randomId(), workflowSources = workflowSources)
  val supervisorProbe: TestProbe = TestProbe("supervisorProbe")
  val deathwatch: TestProbe = TestProbe("deathwatch")
  val initializationProbe: TestProbe = TestProbe("initializationProbe")
  val executionProbe: TestProbe = TestProbe("executionProbe")
  val finalizationProbe: TestProbe = TestProbe("finalizationProbe")
  var copyWorkflowLogsProbe: TestProbe = _
  val workflowCallbackProbe: TestProbe = TestProbe("workflowCallbackProbe")
  val AwaitAlmostNothing: FiniteDuration = 100.milliseconds
  val initialJobCtByRootWf = new AtomicInteger()
  val callCachingEnabled = true
  val invalidateBadCacheResults = true
  val groupMetricsProbe: TestProbe = TestProbe("groupMetricsProbe")

  before {
    currentWorkflowId = WorkflowId.randomId()
    copyWorkflowLogsProbe = TestProbe(s"copyWorkflowLogsProbe-$currentWorkflowId")
    // Clear the supervisor probe of anything remaining from previous runs:
    supervisorProbe.receiveWhile(max = 1.second, idle = 1.second) { case _ =>
      println("Ignoring excess message to WMA: ")
    }
  }

  private val workflowHeartbeatConfig = WorkflowHeartbeatConfig(ConfigFactory.load())

  private def createWorkflowActor(state: WorkflowActorState,
                                  extraPathBuilderFactory: Option[PathBuilderFactory] = None,
                                  workflowCallbackActor: Option[ActorRef] = None,
                                  initializationMaxRetries: Int = 3,
                                  initializationInterval: FiniteDuration = 10.millis
  ) = {
    val actor = TestFSMRef(
      factory = new WorkflowActorWithTestAddons(
        finalizationProbe = finalizationProbe,
        workflowId = currentWorkflowId,
        startState = Submitted,
        workflowSources = workflowSources,
        callCachingEnabled = callCachingEnabled,
        invalidateBadCacheResults = invalidateBadCacheResults,
        conf = ConfigFactory.load,
        ioActor = system.actorOf(SimpleIoActor.props, s"ioActor-$currentWorkflowId"),
        serviceRegistryActor = mockServiceRegistryActor,
        workflowLogCopyRouter = copyWorkflowLogsProbe.ref,
        workflowCallbackActor = workflowCallbackActor,
        jobStoreActor = system.actorOf(AlwaysHappyJobStoreActor.props, s"jobStoreActor-$currentWorkflowId"),
        subWorkflowStoreActor =
          system.actorOf(AlwaysHappySubWorkflowStoreActor.props, s"subWorkflowStoreActor-$currentWorkflowId"),
        callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props, s"callCacheReadActor-$currentWorkflowId"),
        callCacheWriteActor = system.actorOf(EmptyCallCacheWriteActor.props, s"callCacheWriteActor-$currentWorkflowId"),
        dockerHashActor = system.actorOf(EmptyDockerHashActor.props, s"dockerHashActor-$currentWorkflowId"),
        jobRestartCheckTokenDispenserActor = TestProbe(s"jobRestartCheckTokenDispenserActor-$currentWorkflowId").ref,
        jobExecutionTokenDispenserActor = TestProbe(s"jobExecutionTokenDispenserActor-$currentWorkflowId").ref,
        workflowStoreActor = system.actorOf(Props.empty, s"workflowStoreActor-$currentWorkflowId"),
        workflowHeartbeatConfig = workflowHeartbeatConfig,
        totalJobsByRootWf = initialJobCtByRootWf,
        extraPathBuilderFactory = extraPathBuilderFactory,
        initializationMaxRetries = initializationMaxRetries,
        initializationInterval = initializationInterval,
        workflowInitializationActorProbe = initializationProbe,
        workflowExecutionActorProbe = executionProbe,
        groupMetricsActor = groupMetricsProbe.ref
      ),
      supervisor = supervisorProbe.ref,
      name = s"workflowActor-$currentWorkflowId"
    )
    actor.setState(
      stateName = state,
      stateData = WorkflowActorData(Option(currentLifecycleActor.ref),
                                    Option(descriptor),
                                    AllBackendInitializationData.empty,
                                    StateCheckpoint(InitializingWorkflowState),
                                    Submitted
      )
    )
    actor
  }

  implicit val TimeoutDuration: FiniteDuration = CromwellTestKitSpec.TimeoutDuration

  private def workflowManagerActorExpectsSingleWorkCompleteNotification(endState: WorkflowState) = {
    supervisorProbe.expectMsgPF(TimeoutDuration) {
      case wawc: WorkflowActorWorkComplete => wawc.finalState should be(endState)
      case other => fail(s"Unexpected message to WMA while waiting for work complete: $other")
    }
    supervisorProbe.expectNoMessage(AwaitAlmostNothing)
  }

  "WorkflowActor" should {

    "run Finalization actor if Initialization repeatedly fails" in {
      val actor = createWorkflowActor(InitializingWorkflowState, initializationMaxRetries = 3)
      deathwatch watch actor

      // Expect the WorkflowActor to retry up to three times:
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (1)")))
      initializationProbe.expectMsg(StartInitializationCommand)
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (2)")))
      initializationProbe.expectMsg(StartInitializationCommand)
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (3)")))

      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectMsgPF(TimeoutDuration) { case x: WorkflowFailedResponse =>
        x.workflowId should be(currentWorkflowId)
      }
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowFailed)
      deathwatch.expectTerminated(actor)
    }

    "begin execution if Initialization eventually succeeds" in {
      val actor = createWorkflowActor(InitializingWorkflowState, initializationMaxRetries = 3)
      deathwatch watch actor

      // Expect the WorkflowActor to retry up to three times:
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (1)")))
      initializationProbe.expectMsg(StartInitializationCommand)
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (2)")))
      initializationProbe.expectMsg(StartInitializationCommand)
      actor ! WorkflowInitializationSucceededResponse(AllBackendInitializationData(Map.empty))

      executionProbe.expectMsg(ExecuteWorkflowCommand)
      actor.stateName should be(ExecutingWorkflowState)
      // Tidy up (and satisfy the deathwatch):
      system.stop(actor)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Initialization is aborted" in {
      val actor = createWorkflowActor(InitializingWorkflowState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      eventually(actor.stateName should be(WorkflowAbortingState))
      currentLifecycleActor.expectMsgPF(TimeoutDuration) { case EngineLifecycleActorAbortCommand =>
        actor ! WorkflowInitializationAbortedResponse
      }
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowAborted)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Initialization is aborted during a retry" in {
      val actor = createWorkflowActor(InitializingWorkflowState, initializationInterval = 1.second)
      deathwatch watch actor

      // Set the stage with a few unfortunate retries:
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (1)")))
      initializationProbe.expectMsg(StartInitializationCommand)
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (2)")))
      initializationProbe.expectMsg(StartInitializationCommand)

      actor ! AbortWorkflowCommand
      eventually(actor.stateName should be(WorkflowAbortingState))
      // Because we failed a few times, the actor test's initial "currentLifecycleActor will have been replaced by this
      // new initializationProbe:
      initializationProbe.expectMsgPF(TimeoutDuration) { case EngineLifecycleActorAbortCommand =>
        actor ! WorkflowInitializationAbortedResponse
      }

      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowAborted)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Initialization is aborted between retries" in {
      // Set the interval long so that we can guarantee the "between retries" part
      val actor = createWorkflowActor(InitializingWorkflowState, initializationInterval = 10.seconds)
      deathwatch watch actor

      // Set the stage with a few unfortunate retries:
      actor ! WorkflowInitializationFailedResponse(Seq(new Exception("Initialization Failed (1)")))
      eventually(actor.stateData.currentLifecycleStateActor should be(None))

      actor ! AbortWorkflowCommand
      // Because there are no active lifecycle actors, this actor should jump to Finalizing without
      // needing any further input:
      eventually(actor.stateName should be(FinalizingWorkflowState))

      // Expect the mailboxes for the initialization actor to be empty (and check "currentLifecycleActor" for good measure)
      currentLifecycleActor.expectNoMessage(10.millis)
      initializationProbe.expectNoMessage(10.millis)

      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowAborted)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization if Execution fails" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch actor
      actor ! WorkflowExecutionFailedResponse(Map.empty, new Exception("Execution Failed"))
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectMsgPF(TimeoutDuration) { case x: WorkflowFailedResponse =>
        x.workflowId should be(currentWorkflowId)
      }
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowFailed)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Execution is aborted" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      eventually(actor.stateName should be(WorkflowAbortingState))
      currentLifecycleActor.expectMsgPF(CromwellTestKitSpec.TimeoutDuration) { case EngineLifecycleActorAbortCommand =>
        actor ! WorkflowExecutionAbortedResponse(Map.empty)
      }
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowAborted)
      deathwatch.expectTerminated(actor)
    }

    "run Finalization actor if Execution succeeds" in {
      val actor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch actor
      actor ! WorkflowExecutionSucceededResponse(Map.empty, Set(WorkflowId.randomId()), CallOutputs.empty)
      finalizationProbe.expectMsg(StartFinalizationCommand)
      actor.stateName should be(FinalizingWorkflowState)
      actor ! WorkflowFinalizationSucceededResponse
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowSucceeded)
      deathwatch.expectTerminated(actor)
    }

    "not run Finalization actor if aborted when in WorkflowUnstartedState" in {
      val actor = createWorkflowActor(WorkflowUnstartedState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowAborted)
      deathwatch.expectTerminated(actor)
    }

    "not run Finalization actor if aborted when in MaterializingWorkflowDescriptorState" in {
      val actor = createWorkflowActor(MaterializingWorkflowDescriptorState)
      deathwatch watch actor
      actor ! AbortWorkflowCommand
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowAborted)
      deathwatch.expectTerminated(actor)
    }

    "copy workflow logs in the event of MaterializeWorkflowDescriptorFailureResponse" in {
      val actor = createWorkflowActor(MaterializingWorkflowDescriptorState)
      deathwatch watch actor

      copyWorkflowLogsProbe.expectNoMessage(AwaitAlmostNothing)
      actor ! MaterializeWorkflowDescriptorFailureResponse(
        new Exception("Intentionally failing workflow materialization to test log copying")
      )
      copyWorkflowLogsProbe.expectMsg(CopyWorkflowLogsActor.Copy(currentWorkflowId, mockDir))
      supervisorProbe.expectMsgPF(TimeoutDuration) { case _: WorkflowFailedResponse => /* success! */ }
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowFailed)
      deathwatch.expectTerminated(actor)
    }

    "abort execution before escalating failure if one of its child actors crashed" in {
      val workflowActor = createWorkflowActor(ExecutingWorkflowState)
      deathwatch watch workflowActor

      workflowActor.children.head ! Kill

      eventually(workflowActor.stateName should be(WorkflowAbortingState))
      currentLifecycleActor.expectMsgPF(TimeoutDuration) { case EngineLifecycleActorAbortCommand =>
        workflowActor ! WorkflowExecutionAbortedResponse(Map.empty)
      }
      finalizationProbe.expectMsg(StartFinalizationCommand)
      workflowActor.stateName should be(FinalizingWorkflowState)
      workflowActor ! WorkflowFinalizationSucceededResponse
      supervisorProbe.expectMsgPF(TimeoutDuration) { case _: WorkflowFailedResponse => /* success! */ }
      workflowManagerActorExpectsSingleWorkCompleteNotification(WorkflowFailed)
      deathwatch.expectTerminated(workflowActor)
    }

    "log an error when a path builder factory initialization fails" in {
      EventFilter
        .error(start = "Failed to copy workflow log", pattern = ".*Failing as requested.*", occurrences = 1)
        .intercept {
          val _ = createWorkflowActor(WorkflowSucceededState, Option(new FailingPathBuilderFactory()))
        }
    }

    "log an error when a path builder factory initialization throws" in {
      EventFilter
        .error(start = "Failed to copy workflow log", pattern = ".*Throwing as requested.*", occurrences = 1)
        .intercept {
          val _ = createWorkflowActor(WorkflowSucceededState, Option(new ThrowingPathBuilderFactory()))
        }
    }

    "send a workflow callback message" in {
      val actor = createWorkflowActor(ExecutingWorkflowState, workflowCallbackActor = Option(workflowCallbackProbe.ref))
      deathwatch watch actor
      val mockOutputs = WomMocks.mockOutputExpectations(Map("foo" -> WomString("bar")))
      val msg = WorkflowCallbackActor.PerformCallbackCommand(currentWorkflowId,
                                                             Some(mockUri),
                                                             WorkflowSucceeded,
                                                             mockOutputs,
                                                             List.empty
      )

      workflowCallbackProbe.expectNoMessage(AwaitAlmostNothing)
      actor ! WorkflowExecutionSucceededResponse(Map.empty, Set(currentWorkflowId), mockOutputs, Set.empty)
      actor ! WorkflowFinalizationSucceededResponse
      workflowCallbackProbe.expectMsg(msg)
      deathwatch.expectTerminated(actor)
    }

    "send a workflow callback message for a failing workflow" in {
      val actor =
        createWorkflowActor(FinalizingWorkflowState, workflowCallbackActor = Option(workflowCallbackProbe.ref))
      deathwatch watch actor
      val errorText = "oh nooo :("
      val msg = WorkflowCallbackActor.PerformCallbackCommand(currentWorkflowId,
                                                             Some(mockUri),
                                                             WorkflowFailed,
                                                             CallOutputs.empty,
                                                             List(errorText)
      )

      workflowCallbackProbe.expectNoMessage(AwaitAlmostNothing)
      actor ! WorkflowFinalizationFailedResponse(Seq(new RuntimeException(errorText)))
      workflowCallbackProbe.expectMsg(msg)
      deathwatch.expectTerminated(actor)
    }
  }
}

class FailingPathBuilderFactory() extends PathBuilderFactory {
  override def withOptions(
    options: WorkflowOptions
  )(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] =
    Future(throw new Exception("Failing as requested"))
}

class ThrowingPathBuilderFactory() extends PathBuilderFactory {
  override def withOptions(
    options: WorkflowOptions
  )(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] =
    throw new Exception("Throwing as requested")
}

class WorkflowActorWithTestAddons(val finalizationProbe: TestProbe,
                                  workflowId: WorkflowId,
                                  startState: StartableState,
                                  workflowSources: WorkflowSourceFilesCollection,
                                  conf: Config,
                                  callCachingEnabled: Boolean,
                                  invalidateBadCacheResults: Boolean,
                                  ioActor: ActorRef,
                                  serviceRegistryActor: ActorRef,
                                  workflowLogCopyRouter: ActorRef,
                                  workflowCallbackActor: Option[ActorRef],
                                  jobStoreActor: ActorRef,
                                  subWorkflowStoreActor: ActorRef,
                                  callCacheReadActor: ActorRef,
                                  callCacheWriteActor: ActorRef,
                                  dockerHashActor: ActorRef,
                                  jobRestartCheckTokenDispenserActor: ActorRef,
                                  jobExecutionTokenDispenserActor: ActorRef,
                                  workflowStoreActor: ActorRef,
                                  workflowHeartbeatConfig: WorkflowHeartbeatConfig,
                                  totalJobsByRootWf: AtomicInteger,
                                  extraPathBuilderFactory: Option[PathBuilderFactory],
                                  initializationMaxRetries: Int,
                                  initializationInterval: FiniteDuration,
                                  workflowInitializationActorProbe: TestProbe,
                                  workflowExecutionActorProbe: TestProbe,
                                  groupMetricsActor: ActorRef
) extends WorkflowActor(
      workflowToStart = WorkflowToStart(id = workflowId,
                                        submissionTime = OffsetDateTime.now,
                                        state = startState,
                                        sources = workflowSources,
                                        hogGroup = HogGroup("foo")
      ),
      conf = conf,
      callCachingEnabled = callCachingEnabled,
      invalidateBadCacheResults = invalidateBadCacheResults,
      ioActor = ioActor,
      serviceRegistryActor = serviceRegistryActor,
      workflowLogCopyRouter = workflowLogCopyRouter,
      workflowCallbackActor,
      jobStoreActor = jobStoreActor,
      subWorkflowStoreActor = subWorkflowStoreActor,
      callCacheReadActor = callCacheReadActor,
      callCacheWriteActor = callCacheWriteActor,
      dockerHashActor = dockerHashActor,
      jobRestartCheckTokenDispenserActor = jobRestartCheckTokenDispenserActor,
      jobExecutionTokenDispenserActor = jobExecutionTokenDispenserActor,
      backendSingletonCollection = BackendSingletonCollection(Map.empty),
      workflowStoreActor = workflowStoreActor,
      serverMode = true,
      workflowHeartbeatConfig = workflowHeartbeatConfig,
      totalJobsByRootWf = totalJobsByRootWf,
      fileHashCacheActorProps = None,
      blacklistCache = None,
      groupMetricsActor = groupMetricsActor
    ) {

  override def createInitializationActor(workflowDescriptor: EngineWorkflowDescriptor, name: String): ActorRef =
    workflowInitializationActorProbe.ref
  override def createWorkflowExecutionActor(workflowDescriptor: EngineWorkflowDescriptor,
                                            data: WorkflowActorData
  ): ActorRef = workflowExecutionActorProbe.ref

  override val initializationRetryInterval: FiniteDuration = initializationInterval
  override val maxInitializationAttempts: Int = initializationMaxRetries

  override val pathBuilderFactories: List[PathBuilderFactory] = extraPathBuilderFactory match {
    case Some(pbf) => EngineFilesystems.configuredPathBuilderFactories :+ pbf
    case None => EngineFilesystems.configuredPathBuilderFactories
  }

  override def makeFinalizationActor(workflowDescriptor: EngineWorkflowDescriptor,
                                     jobExecutionMap: JobExecutionMap,
                                     workflowOutputs: CallOutputs
  ): ActorRef = finalizationProbe.ref
}
