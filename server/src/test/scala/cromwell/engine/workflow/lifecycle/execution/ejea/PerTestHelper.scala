package cromwell.engine.workflow.lifecycle.execution.ejea

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend.BackendJobExecutionActor.{ExecuteJobCommand, RecoverJobCommand}
import cromwell.backend._
import cromwell.backend.standard.callcaching._
import cromwell.core.callcaching._
import cromwell.core.{CallOutputs, HogGroup, WorkflowId, WorkflowOptions}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor.{
  CallCachingParameters,
  EJEAData,
  EngineJobExecutionActorState,
  ResponsePendingData
}
import cromwell.engine.workflow.mocks.{DeclarationMock, TaskMock, WdlWomExpressionMock}
import cromwell.services.CallCaching.CallCachingEntryId
import cromwell.util.AkkaTestUtil._
import cromwell.util.WomMocks
import wom.callable.Callable.{OutputDefinition, OverridableInputDefinitionWithDefault}
import wom.callable.CallableTaskDefinition
import wom.expression.{IoFunctionSet, NoIoFunctionSet}
import wom.graph.{CommandCallNode, WomIdentifier}
import wom.types.{WomIntegerType, WomStringType}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration.FiniteDuration
import scala.util.{Success, Try}

private[ejea] class PerTestHelper(implicit val system: ActorSystem)
    extends TaskMock
    with WdlWomExpressionMock
    with DeclarationMock {

  val workflowId: WorkflowId = WorkflowId.randomId()
  val workflowName = "wf"
  val taskName = "foobar"
  val jobFqn = s"$workflowName.$taskName"
  val jobIndex: Option[Int] = Option(1)
  val jobAttempt = 1

  val task: CallableTaskDefinition = WomMocks
    .mockTaskDefinition(taskName)
    .copy(
      inputs = List(OverridableInputDefinitionWithDefault("inInt", WomIntegerType, mockIntExpression(543))),
      outputs = List(OutputDefinition("outString", WomStringType, mockStringExpression("hello")))
    )

  val call: CommandCallNode = WomMocks.mockTaskCall(WomIdentifier(taskName, jobFqn), task)
  val jobDescriptorKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, jobIndex, jobAttempt)

  val backendWorkflowDescriptor: BackendWorkflowDescriptor = BackendWorkflowDescriptor(id = workflowId,
                                                                                       callable = null,
                                                                                       knownValues = null,
                                                                                       workflowOptions =
                                                                                         WorkflowOptions.empty,
                                                                                       customLabels = null,
                                                                                       hogGroup = HogGroup("foo"),
                                                                                       List.empty,
                                                                                       None
  )
  val backendJobDescriptor: BackendJobDescriptor =
    BackendJobDescriptor(
      workflowDescriptor = backendWorkflowDescriptor,
      key = jobDescriptorKey,
      runtimeAttributes = Map.empty,
      evaluatedTaskInputs = Map.empty,
      maybeCallCachingEligible = FloatingDockerTagWithoutHash("ubuntu:latest"),
      dockerSize = None,
      prefetchedKvStoreEntries = Map.empty
    )

  var fetchCachedResultsActorCreations: ExpectOne[(CallCachingEntryId, Seq[OutputDefinition])] = NothingYet
  var jobHashingInitializations: ExpectOne[(BackendJobDescriptor, CallCachingActivity)] = NothingYet
  var invalidateCacheActorCreations: ExpectOne[CallCachingEntryId] = NothingYet

  val deathwatch: TestProbe = TestProbe()
  val bjeaProbe: TestProbe = TestProbe()
  val bjeaProps: Props = bjeaProbe.props
  val replyToProbe: TestProbe = TestProbe()
  val parentProbe: TestProbe = TestProbe()
  val serviceRegistryProbe: TestProbe = TestProbe()
  val ioActorProbe: TestProbe = TestProbe()
  val jobStoreProbe: TestProbe = TestProbe()
  val callCacheReadActorProbe: TestProbe = TestProbe()
  val callCacheWriteActorProbe: TestProbe = TestProbe()
  val dockerHashActorProbe: TestProbe = TestProbe()
  val callCacheHitCopyingProbe: TestProbe = TestProbe()
  val jobPreparationProbe: TestProbe = TestProbe()
  val jobRestartCheckTokenDispenserProbe: TestProbe = TestProbe()
  val jobExecutionTokenDispenserProbe: TestProbe = TestProbe()
  val ejhaProbe: TestProbe = TestProbe()
  val groupMetricsProbe: TestProbe = TestProbe()

  def buildFactory(backendConfigurationDescriptor: BackendConfigurationDescriptor): BackendLifecycleActorFactory =
    new BackendLifecycleActorFactory {

      override val name = "PerTestHelper"

      override val configurationDescriptor: BackendConfigurationDescriptor = backendConfigurationDescriptor

      override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                          initializationData: Option[BackendInitializationData],
                                          serviceRegistryActor: ActorRef,
                                          ioActor: ActorRef,
                                          backendSingletonActor: Option[ActorRef],
                                          groupMetricsActor: ActorRef
      ): Props = bjeaProps

      override def cacheHitCopyingActorProps: Option[
        (BackendJobDescriptor,
         Option[BackendInitializationData],
         ActorRef,
         ActorRef,
         Int,
         Option[BlacklistCache]
        ) => Props
      ] = Option((_, _, _, _, _, _) => callCacheHitCopyingProbe.props)

      override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                               jobKey: BackendJobDescriptorKey,
                                               initializationData: Option[BackendInitializationData],
                                               ioActorProxy: ActorRef,
                                               ec: ExecutionContext
      ): IoFunctionSet =
        NoIoFunctionSet

      override def fileHashingActorProps: Option[
        (BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef, Option[ActorRef]) => Props
      ] =
        Option((_, _, _, _, _) => Props.empty)

      // These two factory methods should never be called from EJEA or any of its descendants:
      override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                  ioActor: ActorRef,
                                                  calls: Set[CommandCallNode],
                                                  jobExecutionMap: JobExecutionMap,
                                                  workflowOutputs: CallOutputs,
                                                  initializationData: Option[BackendInitializationData]
      ): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
      override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                    ioActor: ActorRef,
                                                    calls: Set[CommandCallNode],
                                                    serviceRegistryActor: ActorRef,
                                                    restarting: Boolean
      ): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
    }

  def buildEJEA(restarting: Boolean = true,
                callCachingMode: CallCachingMode = CallCachingOff,
                callCachingMaxFailedCopyAttempts: Int = 1000000,
                backendConfigurationDescriptor: BackendConfigurationDescriptor = TestConfig.emptyBackendConfigDescriptor
  )(implicit
    startingState: EngineJobExecutionActorState
  ): TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea] = {

    val factory: BackendLifecycleActorFactory = buildFactory(backendConfigurationDescriptor)
    val descriptor = EngineWorkflowDescriptor(WomMocks.mockWorkflowDefinition(workflowName),
                                              backendWorkflowDescriptor,
                                              null,
                                              null,
                                              null,
                                              callCachingMode
    )

    val callCachingParameters = CallCachingParameters(
      mode = callCachingMode,
      readActor = callCacheReadActorProbe.ref,
      writeActor = callCacheWriteActorProbe.ref,
      fileHashCacheActor = Option(dockerHashActorProbe.ref),
      maxFailedCopyAttempts = callCachingMaxFailedCopyAttempts,
      blacklistCache = None,
      fileHashBatchSize = 100
    )

    val myBrandNewEjea = new TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea](
      system,
      Props(
        new MockEjea(
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
          dockerHashActor = dockerHashActorProbe.ref,
          jobRestartCheckTokenDispenserActor = jobRestartCheckTokenDispenserProbe.ref,
          jobExecutionTokenDispenserActor = jobExecutionTokenDispenserProbe.ref,
          callCachingParameters = callCachingParameters,
          groupMetricsActor = groupMetricsProbe.ref
        )
      ),
      parentProbe.ref,
      s"EngineJobExecutionActorSpec-$workflowId"
    )

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
                             dockerHashActor: ActorRef,
                             jobRestartCheckTokenDispenserActor: ActorRef,
                             jobExecutionTokenDispenserActor: ActorRef,
                             callCachingParameters: EngineJobExecutionActor.CallCachingParameters,
                             groupMetricsActor: ActorRef
) extends EngineJobExecutionActor(
      replyTo = replyTo,
      jobDescriptorKey = jobDescriptorKey,
      workflowDescriptor = workflowDescriptor,
      backendLifecycleActorFactory = factory,
      initializationData = initializationData,
      restarting = restarting,
      serviceRegistryActor = serviceRegistryActor,
      ioActor = ioActor,
      jobStoreActor = jobStoreActor,
      workflowDockerLookupActor = dockerHashActor,
      jobRestartCheckTokenDispenserActor = jobRestartCheckTokenDispenserActor,
      jobExecutionTokenDispenserActor = jobExecutionTokenDispenserActor,
      backendSingletonActor = None,
      command = if (restarting) RecoverJobCommand else ExecuteJobCommand,
      callCachingParameters = callCachingParameters,
      groupMetricsActor = groupMetricsActor
    ) {

  implicit val system: ActorSystem = context.system
  override def makeFetchCachedResultsActor(cacheId: CallCachingEntryId): Unit =
    helper.fetchCachedResultsActorCreations = helper.fetchCachedResultsActorCreations.foundOne((cacheId, null))
  override def initializeJobHashing(jobDescriptor: BackendJobDescriptor,
                                    activity: CallCachingActivity,
                                    callCachingEligible: CallCachingEligible
  ): Try[ActorRef] = {
    helper.jobHashingInitializations = helper.jobHashingInitializations.foundOne((jobDescriptor, activity))
    Success(helper.ejhaProbe.ref)
  }
  override def createBackendJobExecutionActor(data: ResponsePendingData): ActorRef = helper.bjeaProbe.ref
  override def invalidateCacheHit(cacheId: CallCachingEntryId): Unit = helper.invalidateCacheActorCreations =
    helper.invalidateCacheActorCreations.foundOne(cacheId)
  override def createJobPreparationActor(jobPrepProps: Props, name: String): ActorRef = jobPreparationProbe.ref
  override def onTimedTransition(from: EngineJobExecutionActorState,
                                 to: EngineJobExecutionActorState,
                                 duration: FiniteDuration
  ): Unit = {}
}
