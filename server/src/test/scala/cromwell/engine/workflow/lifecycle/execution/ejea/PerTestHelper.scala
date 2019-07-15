package cromwell.engine.workflow.lifecycle.execution.ejea

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend.BackendJobExecutionActor.{ExecuteJobCommand, RecoverJobCommand}
import cromwell.backend._
import cromwell.backend.standard.callcaching._
import cromwell.core.callcaching._
import cromwell.core.{CallOutputs, HogGroup, WorkflowId, WorkflowOptions}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor.{EJEAData, EngineJobExecutionActorState, ResponsePendingData}
import cromwell.engine.workflow.mocks.{DeclarationMock, TaskMock, WdlWomExpressionMock}
import cromwell.util.AkkaTestUtil._
import cromwell.util.WomMocks
import org.specs2.mock.Mockito
import wom.callable.Callable.{OutputDefinition, OverridableInputDefinitionWithDefault}
import wom.expression.{IoFunctionSet, NoIoFunctionSet}
import wom.graph.{CommandCallNode, WomIdentifier}
import wom.types.{WomIntegerType, WomStringType}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration.FiniteDuration
import scala.util.Success

private[ejea] class PerTestHelper(implicit val system: ActorSystem) extends Mockito with TaskMock with WdlWomExpressionMock with DeclarationMock {

  val workflowId = WorkflowId.randomId()
  val workflowName = "wf"
  val taskName = "foobar"
  val jobFqn = s"$workflowName.$taskName"
  val jobIndex = Some(1)
  val jobAttempt = 1

  val task = WomMocks.mockTaskDefinition(taskName).copy(
    inputs = List(OverridableInputDefinitionWithDefault("inInt", WomIntegerType, mockIntExpression(543))),
    outputs = List(OutputDefinition("outString", WomStringType, mockStringExpression("hello")))
  )

  val call: CommandCallNode = WomMocks.mockTaskCall(WomIdentifier(taskName, jobFqn), task)
  val jobDescriptorKey = BackendJobDescriptorKey(call, jobIndex, jobAttempt)

  val backendWorkflowDescriptor = BackendWorkflowDescriptor(id = workflowId,
    callable = null,
    knownValues = null,
    workflowOptions = WorkflowOptions.empty,
    customLabels = null,
    hogGroup = HogGroup("foo"),
    List.empty,
    None)
  val backendJobDescriptor = BackendJobDescriptor(backendWorkflowDescriptor, jobDescriptorKey, runtimeAttributes = Map.empty, evaluatedTaskInputs = Map.empty, FloatingDockerTagWithoutHash("ubuntu:latest"), None, Map.empty)

  var fetchCachedResultsActorCreations: ExpectOne[(CallCachingEntryId, Seq[OutputDefinition])] = NothingYet
  var jobHashingInitializations: ExpectOne[(BackendJobDescriptor, CallCachingActivity)] = NothingYet
  var invalidateCacheActorCreations: ExpectOne[CallCachingEntryId] = NothingYet

  val deathwatch = TestProbe()
  val bjeaProbe = TestProbe()
  val bjeaProps = bjeaProbe.props
  val replyToProbe = TestProbe()
  val parentProbe = TestProbe()
  val serviceRegistryProbe = TestProbe()
  val ioActorProbe = TestProbe()
  val jobStoreProbe = TestProbe()
  val callCacheReadActorProbe = TestProbe()
  val callCacheWriteActorProbe = TestProbe()
  val dockerHashActorProbe = TestProbe()
  val callCacheHitCopyingProbe = TestProbe()
  val jobPreparationProbe = TestProbe()
  val jobTokenDispenserProbe = TestProbe()
  val ejhaProbe = TestProbe()

  def buildFactory() = new BackendLifecycleActorFactory {

    override val name = "PerTestHelper"

    override val configurationDescriptor = TestConfig.emptyBackendConfigDescriptor

    override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                        initializationData: Option[BackendInitializationData],
                                        serviceRegistryActor: ActorRef,
                                        ioActor: ActorRef,
                                        backendSingletonActor: Option[ActorRef]): Props = bjeaProps

    override def cacheHitCopyingActorProps: Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef, Int, Option[BlacklistCache]) => Props] = Option((_, _, _, _, _, _) => callCacheHitCopyingProbe.props)

    override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                             jobKey: BackendJobDescriptorKey,
                                             initializationData: Option[BackendInitializationData],
                                             ioActorProxy: ActorRef,
                                             ec: ExecutionContext): IoFunctionSet = {
      NoIoFunctionSet
    }

    override def fileHashingActorProps:
    Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef, Option[ActorRef]) => Props] = {
      Option(fileHashingActorInner(classOf[DefaultStandardFileHashingActor]))
    }

    def fileHashingActorInner(standardFileHashingActor: Class[_ <: StandardFileHashingActor])
                             (jobDescriptor: BackendJobDescriptor,
                              initializationDataOption: Option[BackendInitializationData],
                              serviceRegistryActor: ActorRef,
                              ioActor: ActorRef,
                              fileHashCacheActor: Option[ActorRef]): Props = {
      Props.empty
    }

    // These two factory methods should never be called from EJEA or any of its descendants:
    override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                ioActor: ActorRef,
                                                calls: Set[CommandCallNode],
                                                jobExecutionMap: JobExecutionMap,
                                                workflowOutputs: CallOutputs,
                                                initializationData: Option[BackendInitializationData]): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
    override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                  ioActor: ActorRef,
                                                  calls: Set[CommandCallNode],
                                                  serviceRegistryActor: ActorRef,
                                                  restarting: Boolean): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
  }

  def buildEJEA(restarting: Boolean = true,
                callCachingMode: CallCachingMode = CallCachingOff)
               (implicit startingState: EngineJobExecutionActorState): TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea] = {

    val factory: BackendLifecycleActorFactory = buildFactory()
    val descriptor = EngineWorkflowDescriptor(WomMocks.mockWorkflowDefinition(workflowName), backendWorkflowDescriptor, null, null, null, callCachingMode)

    val myBrandNewEjea = new TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea](system, Props(new MockEjea(
      helper = this,
      jobPreparationProbe = jobPreparationProbe,
      replyTo = replyToProbe.ref,
      jobDescriptorKey = jobDescriptorKey,
      workflowDescriptor = descriptor,
      factory = factory,
      initializationData = None,
      restarting = restarting,
      serviceRegistryActor = serviceRegistryProbe.ref,
      ioActor = ioActorProbe.ref,
      jobStoreActor = jobStoreProbe.ref,
      callCacheReadActor = callCacheReadActorProbe.ref,
      callCacheWriteActor = callCacheWriteActorProbe.ref,
      dockerHashActor = dockerHashActorProbe.ref,
      jobTokenDispenserActor = jobTokenDispenserProbe.ref,
      callCachingMode = callCachingMode
    )), parentProbe.ref, s"EngineJobExecutionActorSpec-$workflowId")

    deathwatch watch myBrandNewEjea
    myBrandNewEjea.setStateInline(state = startingState)
  }
}

private[ejea] class MockEjea(helper: PerTestHelper,
                             jobPreparationProbe: TestProbe,
                             replyTo: ActorRef,
                             jobDescriptorKey: BackendJobDescriptorKey,
                             workflowDescriptor: EngineWorkflowDescriptor,
                             factory: BackendLifecycleActorFactory,
                             initializationData: Option[BackendInitializationData],
                             restarting: Boolean,
                             serviceRegistryActor: ActorRef,
                             ioActor: ActorRef,
                             jobStoreActor: ActorRef,
                             callCacheReadActor: ActorRef,
                             callCacheWriteActor: ActorRef,
                             dockerHashActor: ActorRef,
                             jobTokenDispenserActor: ActorRef,
                             callCachingMode: CallCachingMode) extends EngineJobExecutionActor(replyTo, jobDescriptorKey, workflowDescriptor, factory,
  initializationData, restarting, serviceRegistryActor, ioActor,
  jobStoreActor, callCacheReadActor, callCacheWriteActor,
  dockerHashActor, jobTokenDispenserActor, None, callCachingMode,
  if (restarting) RecoverJobCommand else ExecuteJobCommand, fileHashCachingActor = None, blacklistCache = None
) {

  implicit val system = context.system
  override def makeFetchCachedResultsActor(cacheId: CallCachingEntryId) = helper.fetchCachedResultsActorCreations = helper.fetchCachedResultsActorCreations.foundOne((cacheId, null))
  override def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity, callCachingEligible: CallCachingEligible) = {
    helper.jobHashingInitializations = helper.jobHashingInitializations.foundOne((jobDescriptor, activity))
    Success(helper.ejhaProbe.ref)
  }
  override def createBackendJobExecutionActor(data: ResponsePendingData) = helper.bjeaProbe.ref
  override def invalidateCacheHit(cacheId: CallCachingEntryId): Unit = { helper.invalidateCacheActorCreations = helper.invalidateCacheActorCreations.foundOne(cacheId) }
  override def createJobPreparationActor(jobPrepProps: Props, name: String) = jobPreparationProbe.ref
  override def onTimedTransition(from: EngineJobExecutionActorState, to: EngineJobExecutionActorState, duration: FiniteDuration) = {}
}
