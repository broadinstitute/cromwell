package cromwell.engine.workflow.lifecycle.execution.ejea

import java.util.UUID

import _root_.wdl._
import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend._
import cromwell.backend.standard.callcaching._
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.core.callcaching._
import cromwell.core.{JobExecutionToken, NoIoFunctionSet, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{EJEAData, EngineJobExecutionActorState}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.mocks.{DeclarationMock, TaskMock, WdlWomExpressionMock}
import cromwell.util.AkkaTestUtil._
import cromwell.util.WomMocks
import org.specs2.mock.Mockito
import wdl4s.parser.WdlParser.Ast
import wom.callable.Callable.{InputDefinitionWithDefault, OutputDefinition}
import wom.core.CallOutputs
import wom.expression.IoFunctionSet
import wom.graph.{TaskCallNode, WomIdentifier}
import wom.types.{WdlIntegerType, WdlStringType}

import scala.util.Success


private[ejea] class PerTestHelper(implicit val system: ActorSystem) extends Mockito with TaskMock with WdlWomExpressionMock with DeclarationMock {

  val workflowId = WorkflowId.randomId()
  val workflowName = "wf"
  val taskName = "foobar"
  val jobFqn = s"$workflowName.$taskName"
  val jobIndex = Some(1)
  val jobAttempt = 1

  val executionToken = JobExecutionToken(JobExecutionTokenType("test", None), UUID.randomUUID())

  val task = WomMocks.mockTaskDefinition(taskName).copy(
    inputs = List(InputDefinitionWithDefault("inInt", WdlIntegerType, mockIntExpression(543))),
    outputs = List(OutputDefinition("outString", WdlStringType, mockStringExpression("hello")))
  )

  val workflow = new WdlWorkflow(
    unqualifiedName = workflowName,
    workflowOutputWildcards = Seq.empty,
    wdlSyntaxErrorFormatter = mock[WdlSyntaxErrorFormatter],
    meta = Map.empty,
    parameterMeta = Map.empty,
    ast = mock[Ast])
  val call: TaskCallNode = WomMocks.mockTaskCall(WomIdentifier(taskName, jobFqn), task)
  val jobDescriptorKey = BackendJobDescriptorKey(call, jobIndex, jobAttempt)

  val backendWorkflowDescriptor = BackendWorkflowDescriptor(workflowId, null, null, null, null)
  val backendJobDescriptor = BackendJobDescriptor(backendWorkflowDescriptor, jobDescriptorKey, runtimeAttributes = Map.empty, inputDeclarations = Map.empty, FloatingDockerTagWithoutHash("ubuntu:latest"), Map.empty)

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

    override def cacheHitCopyingActorProps: Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef) => Props] = Option((_, _, _, _) => callCacheHitCopyingProbe.props)

    override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor, jobKey: BackendJobDescriptorKey, initializationData: Option[BackendInitializationData]): IoFunctionSet = {
      NoIoFunctionSet
    }

    override def fileHashingActorProps:
    Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef) => Props] = {
      Option(fileHashingActorInner(classOf[DefaultStandardFileHashingActor]))
    }

    def fileHashingActorInner(standardFileHashingActor: Class[_ <: StandardFileHashingActor])
                             (jobDescriptor: BackendJobDescriptor,
                              initializationDataOption: Option[BackendInitializationData],
                              serviceRegistryActor: ActorRef,
                              ioActor: ActorRef): Props = {
      Props.empty
    }

    // These two factory methods should never be called from EJEA or any of its descendants:
    override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                ioActor: ActorRef,
                                                calls: Set[TaskCallNode],
                                                jobExecutionMap: JobExecutionMap,
                                                workflowOutputs: CallOutputs,
                                                initializationData: Option[BackendInitializationData]): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
    override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                  ioActor: ActorRef,
                                                  calls: Set[TaskCallNode],
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
                             backendName: String,
                             callCachingMode: CallCachingMode) extends EngineJobExecutionActor(replyTo, jobDescriptorKey, workflowDescriptor, factory, initializationData, restarting, serviceRegistryActor, ioActor, jobStoreActor, callCacheReadActor, callCacheWriteActor, dockerHashActor, jobTokenDispenserActor, None, backendName, callCachingMode) {

  implicit val system = context.system
  override def makeFetchCachedResultsActor(cacheId: CallCachingEntryId) = helper.fetchCachedResultsActorCreations = helper.fetchCachedResultsActorCreations.foundOne((cacheId, null))
  override def initializeJobHashing(jobDescriptor: BackendJobDescriptor, activity: CallCachingActivity, callCachingEligible: CallCachingEligible) = {
    helper.jobHashingInitializations = helper.jobHashingInitializations.foundOne((jobDescriptor, activity))
    Success(helper.ejhaProbe.ref)
  }
  override def invalidateCacheHit(cacheId: CallCachingEntryId): Unit = { helper.invalidateCacheActorCreations = helper.invalidateCacheActorCreations.foundOne(cacheId) }
  override def createJobPreparationActor(jobPrepProps: Props, name: String) = jobPreparationProbe.ref
}
