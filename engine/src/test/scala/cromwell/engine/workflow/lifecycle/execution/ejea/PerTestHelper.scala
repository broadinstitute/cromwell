package cromwell.engine.workflow.lifecycle.execution.ejea

import java.util.UUID

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.backend._
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.core.callcaching.{CallCachingActivity, CallCachingMode, CallCachingOff}
import cromwell.core.{CallOutputs, JobExecutionToken, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{EJEAData, EngineJobExecutionActorState}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.{EngineJobExecutionActor, WorkflowExecutionActorData}
import cromwell.engine.workflow.mocks.{DeclarationMock, TaskMock, WdlExpressionMock}
import cromwell.util.AkkaTestUtil._
import org.specs2.mock.Mockito
import wdl4s._
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}
import wdl4s.parser.WdlParser.Ast
import wdl4s.types.{WdlIntegerType, WdlStringType}


private[ejea] class PerTestHelper(implicit val system: ActorSystem) extends Mockito with TaskMock with WdlExpressionMock with DeclarationMock {

  val workflowId = WorkflowId.randomId()
  val workflowName = "wf"
  val taskName = "foobar"
  val jobFqn = s"$workflowName.$taskName"
  val jobIndex = Some(1)
  val jobAttempt = 1

  val executionToken = JobExecutionToken(JobExecutionTokenType("test", None), UUID.randomUUID())

  val task = mockTask(
    taskName,
    declarations = Seq(mockDeclaration("inInt", WdlIntegerType, mockIntExpression(543))),
    outputs = Seq(("outString", WdlStringType, mockStringExpression("hello")))
  )

  val workflow = new Workflow(
    unqualifiedName = workflowName,
    workflowOutputWildcards = Seq.empty,
    wdlSyntaxErrorFormatter = mock[WdlSyntaxErrorFormatter],
    meta = Map.empty,
    parameterMeta = Map.empty,
    ast = mock[Ast])
  val call: TaskCall = TaskCall(None, task, Map.empty, mock[Ast])
  call.parent_=(workflow)
  val jobDescriptorKey = BackendJobDescriptorKey(call, jobIndex, jobAttempt)

  val backendWorkflowDescriptor = BackendWorkflowDescriptor(workflowId, null, null, null)
  val backendJobDescriptor = BackendJobDescriptor(backendWorkflowDescriptor, jobDescriptorKey, runtimeAttributes = Map.empty, inputDeclarations = Map.empty)

  var fetchCachedResultsActorCreations: ExpectOne[(CallCachingEntryId, Seq[TaskOutput])] = NothingYet
  var jobHashingInitializations: ExpectOne[(BackendJobDescriptor, CallCachingActivity)] = NothingYet
  var callCacheWriteActorCreations: ExpectOne[(CallCacheHashes, JobSucceededResponse)] = NothingYet
  var invalidateCacheActorCreations: ExpectOne[CallCachingEntryId] = NothingYet

  val deathwatch = TestProbe()
  val bjeaProbe = TestProbe()
  val bjeaProps = bjeaProbe.props
  val replyToProbe = TestProbe()
  val parentProbe = TestProbe()
  val serviceRegistryProbe = TestProbe()
  val jobStoreProbe = TestProbe()
  val callCacheReadActorProbe = TestProbe()
  val callCacheHitCopyingProbe = TestProbe()
  val jobPreparationProbe = TestProbe()
  val jobTokenDispenserProbe = TestProbe()

  def buildFactory() = new BackendLifecycleActorFactory {

    override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                        initializationData: Option[BackendInitializationData],
                                        serviceRegistryActor: ActorRef,
                                        backendSingletonActor: Option[ActorRef]): Props = bjeaProps

    override def cacheHitCopyingActorProps: Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef) => Props] = Option((_, _, _) => callCacheHitCopyingProbe.props)

    override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor, jobKey: BackendJobDescriptorKey, initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
      NoFunctions
    }

    // These two factory methods should never be called from EJEA or any of its descendants:
    override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Set[TaskCall],
                                                jobExecutionMap: JobExecutionMap,
                                                workflowOutputs: CallOutputs,
                                                initializationData: Option[BackendInitializationData]): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
    override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                  calls: Set[TaskCall],
                                                  serviceRegistryActor: ActorRef): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
  }

  def buildEJEA(restarting: Boolean = true,
                callCachingMode: CallCachingMode = CallCachingOff)
               (implicit startingState: EngineJobExecutionActorState): TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea] = {

    val factory: BackendLifecycleActorFactory = buildFactory()
    val descriptor = EngineWorkflowDescriptor(mock[WdlNamespaceWithWorkflow], backendWorkflowDescriptor, null, null, null, callCachingMode)

    val myBrandNewEjea = new TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea](system, Props(new MockEjea(
      helper = this,
      jobPreparationProbe = jobPreparationProbe,
      replyTo = replyToProbe.ref,
      jobDescriptorKey = jobDescriptorKey,
      executionData = WorkflowExecutionActorData.empty(descriptor),
      factory = factory,
      initializationData = None,
      restarting = restarting,
      serviceRegistryActor = serviceRegistryProbe.ref,
      jobStoreActor = jobStoreProbe.ref,
      callCacheReadActor = callCacheReadActorProbe.ref,
      jobTokenDispenserActor = jobTokenDispenserProbe.ref,
      backendName = "NOT USED",
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
                             executionData: WorkflowExecutionActorData,
                             factory: BackendLifecycleActorFactory,
                             initializationData: Option[BackendInitializationData],
                             restarting: Boolean,
                             serviceRegistryActor: ActorRef,
                             jobStoreActor: ActorRef,
                             callCacheReadActor: ActorRef,
                             jobTokenDispenserActor: ActorRef,
                             backendName: String,
                             callCachingMode: CallCachingMode) extends EngineJobExecutionActor(replyTo, jobDescriptorKey, executionData, factory, initializationData, restarting, serviceRegistryActor, jobStoreActor, callCacheReadActor, jobTokenDispenserActor, None, backendName, callCachingMode) {

  override def makeFetchCachedResultsActor(cacheId: CallCachingEntryId, taskOutputs: Seq[TaskOutput]) = helper.fetchCachedResultsActorCreations = helper.fetchCachedResultsActorCreations.foundOne((cacheId, taskOutputs))
  override def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity) = helper.jobHashingInitializations = helper.jobHashingInitializations.foundOne((jobDescriptor, activity))
  override def createSaveCacheResultsActor(hashes: CallCacheHashes, success: JobSucceededResponse) = helper.callCacheWriteActorCreations = helper.callCacheWriteActorCreations.foundOne((hashes, success))
  override def invalidateCacheHit(cacheId: CallCachingEntryId): Unit = { helper.invalidateCacheActorCreations = helper.invalidateCacheActorCreations.foundOne(cacheId) }
  override def createJobPreparationActor(jobPrepProps: Props, name: String) = jobPreparationProbe.ref
}
