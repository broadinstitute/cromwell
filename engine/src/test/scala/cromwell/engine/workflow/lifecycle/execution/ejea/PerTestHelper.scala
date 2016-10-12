package cromwell.engine.workflow.lifecycle.execution.ejea

import java.util.UUID

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend.BackendJobExecutionActor.SucceededResponse
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory, BackendWorkflowDescriptor}
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.core.callcaching.{CallCachingActivity, CallCachingMode, CallCachingOff}
import cromwell.core.{ExecutionStore, JobExecutionToken, OutputStore, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.engine.workflow.lifecycle.execution.{EngineJobExecutionActor, WorkflowExecutionActorData}
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{EJEAData, EngineJobExecutionActorState}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes
import org.specs2.mock.Mockito
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s.expression.{NoFunctions, WdlFunctions, WdlStandardLibraryFunctions}
import wdl4s.types.{WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlInteger, WdlString, WdlValue}
import wdl4s._
import cromwell.util.AkkaTestUtil._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._

import scala.util.Success


private[ejea] class PerTestHelper(implicit val system: ActorSystem) extends Mockito {

  val workflowId = WorkflowId.randomId()
  val workflowName = "wf"
  val taskName = "foobar"
  val jobFqn = s"$workflowName.$taskName"
  val jobIndex = Some(1)
  val jobAttempt = 1

  val executionToken = JobExecutionToken(JobExecutionTokenType("test", None), UUID.randomUUID())

  val task = mock[Task]
  task.declarations returns Seq.empty
  task.runtimeAttributes returns RuntimeAttributes(Map.empty)
  task.commandTemplateString returns "!!shazam!!"
  val stringOutputExpression = mock[WdlExpression]
  stringOutputExpression.valueString returns "hello"
  stringOutputExpression.evaluate(any[ScopedLookupFunction], any[ WdlFunctions[WdlValue]]) returns Success(WdlString("hello"))
  task.outputs returns Seq(TaskOutput("outString", WdlStringType, stringOutputExpression))

  val intInputExpression = mock[WdlExpression]
  intInputExpression.valueString returns "543"
  intInputExpression.evaluate(any[ScopedLookupFunction], any[WdlFunctions[WdlValue]]) returns Success(WdlInteger(543))

  val intInputDeclaration = mock[Declaration]
  intInputDeclaration.name returns "inInt"
  intInputDeclaration.expression returns Option(intInputExpression)
  intInputDeclaration.wdlType returns WdlIntegerType
  task.declarations returns Seq(intInputDeclaration)

  val call: Call = Call(None, jobFqn, task, Set.empty, Map.empty, None)
  val jobDescriptorKey = BackendJobDescriptorKey(call, jobIndex, jobAttempt)

  val backendWorkflowDescriptor = BackendWorkflowDescriptor(workflowId, null, null, null)
  val backendJobDescriptor = BackendJobDescriptor(backendWorkflowDescriptor, jobDescriptorKey, runtimeAttributes = Map.empty, inputs = Map.empty)

  var fetchCachedResultsActorCreations: ExpectOne[(CallCachingEntryId, Seq[TaskOutput])] = NothingYet
  var jobHashingInitializations: ExpectOne[(BackendJobDescriptor, CallCachingActivity)] = NothingYet
  var callCacheWriteActorCreations: ExpectOne[(CallCacheHashes, SucceededResponse)] = NothingYet
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
                                                calls: Seq[Call],
                                                executionStore: ExecutionStore,
                                                outputStore: OutputStore,
                                                initializationData: Option[BackendInitializationData]): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
    override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                  calls: Seq[Call],
                                                  serviceRegistryActor: ActorRef): Option[Props] = throw new UnsupportedOperationException("Unexpected finalization actor creation!")
  }

  def buildEJEA(restarting: Boolean = true,
                callCachingMode: CallCachingMode = CallCachingOff)
               (implicit startingState: EngineJobExecutionActorState): TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea] = {

    val factory: BackendLifecycleActorFactory = buildFactory()
    val descriptor = EngineWorkflowDescriptor(backendWorkflowDescriptor, Map.empty, null, null, null, callCachingMode)

    val myBrandNewEjea = new TestFSMRef[EngineJobExecutionActorState, EJEAData, MockEjea](system, Props(new MockEjea(
      helper = this,
      jobPreparationProbe = jobPreparationProbe,
      replyTo = replyToProbe.ref,
      jobDescriptorKey = jobDescriptorKey,
      executionData = WorkflowExecutionActorData(descriptor, ExecutionStore(Map.empty), Map.empty, OutputStore(Map.empty)),
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
  override def createSaveCacheResultsActor(hashes: CallCacheHashes, success: SucceededResponse) = helper.callCacheWriteActorCreations = helper.callCacheWriteActorCreations.foundOne((hashes, success))
  override def invalidateCacheHit(cacheId: CallCachingEntryId): Unit = {
    helper.invalidateCacheActorCreations = helper.invalidateCacheActorCreations.foundOne(cacheId)
  }
  override def createJobPreparationActor(jobPrepProps: Props, name: String) = jobPreparationProbe.ref
}
