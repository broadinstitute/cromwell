package cromwell.engine.workflow.lifecycle.execution.job.preparation

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.TestProbe
import common.mock.MockSugar
import common.validation.ErrorOr.ErrorOr
import cromwell.backend._
import cromwell.core.WorkflowId
import cromwell.docker.DockerMirroring
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData
import cromwell.engine.workflow.lifecycle.execution.job.preparation.JobPreparationTestHelper._
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.services.keyvalue.KeyValueServiceActor.{KvJobKey, ScopedKey}
import wdl.draft2.model.LocallyQualifiedName
import wom.expression.NoIoFunctionSet
import wom.graph.{CommandCallNode, WomIdentifier}
import wom.values.{WomEvaluatedCallInputs, WomValue}

import scala.concurrent.duration.FiniteDuration

class JobPreparationTestHelper(implicit val system: ActorSystem) extends MockSugar {
  val executionData: WorkflowExecutionActorData = mock[WorkflowExecutionActorData]
  val workflowDescriptor: EngineWorkflowDescriptor = mock[EngineWorkflowDescriptor]
  val backendDescriptor: BackendWorkflowDescriptor = mock[BackendWorkflowDescriptor]
  val workflowId: WorkflowId = WorkflowId.randomId()
  workflowDescriptor.backendDescriptor returns backendDescriptor
  workflowDescriptor.id returns workflowId
  workflowDescriptor.possiblyNotRootWorkflowId returns workflowId.toPossiblyNotRoot
  workflowDescriptor.rootWorkflowId returns workflowId.toRoot
  workflowDescriptor.rootWorkflow returns workflowDescriptor
  executionData.workflowDescriptor returns workflowDescriptor
  val call: CommandCallNode =
    CommandCallNode(WomIdentifier("JobPreparationSpec_call"), null, null, null, Set.empty, null, None)
  val mockJobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)
  val serviceRegistryProbe: TestProbe = TestProbe()
  val ioActor: TestProbe = TestProbe()
  val workflowDockerLookupActor: TestProbe = TestProbe()
  val groupMetricsActor: TestProbe = TestProbe()

  val scopedKeyMaker: ScopedKeyMaker = key =>
    ScopedKey(workflowId, KvJobKey("correct.horse.battery.staple", None, 1), key)

  def buildTestJobPreparationActor(
    backpressureTimeout: FiniteDuration,
    noResponseTimeout: FiniteDuration,
    dockerHashCredentials: List[Any],
    inputsAndAttributes: Option[ErrorOr[(WomEvaluatedCallInputs, Map[LocallyQualifiedName, WomValue])]],
    kvStoreKeysForPrefetch: List[String],
    jobKey: BackendJobDescriptorKey = mockJobKey,
    dockerMirroring: Option[DockerMirroring] = None
  ): Props =
    Props(
      new TestJobPreparationActor(
        kvStoreKeysForPrefetch = kvStoreKeysForPrefetch,
        dockerHashCredentialsInput = dockerHashCredentials,
        backpressureWaitTimeInput = backpressureTimeout,
        dockerNoResponseTimeoutInput = noResponseTimeout,
        dockerMirroringInput = dockerMirroring,
        inputsAndAttributes = inputsAndAttributes,
        workflowDescriptor = workflowDescriptor,
        jobKey = jobKey,
        workflowDockerLookupActor = workflowDockerLookupActor.ref,
        serviceRegistryActor = serviceRegistryProbe.ref,
        ioActor = ioActor.ref,
        scopedKeyMaker,
        groupMetricsActor.ref
      )
    )
}

private[preparation] class TestJobPreparationActor(
  kvStoreKeysForPrefetch: List[String],
  dockerHashCredentialsInput: List[Any],
  backpressureWaitTimeInput: FiniteDuration,
  dockerNoResponseTimeoutInput: FiniteDuration,
  dockerMirroringInput: Option[DockerMirroring],
  inputsAndAttributes: Option[ErrorOr[(WomEvaluatedCallInputs, Map[LocallyQualifiedName, WomValue])]],
  workflowDescriptor: EngineWorkflowDescriptor,
  jobKey: BackendJobDescriptorKey,
  workflowDockerLookupActor: ActorRef,
  serviceRegistryActor: ActorRef,
  ioActor: ActorRef,
  scopedKeyMaker: ScopedKeyMaker,
  groupMetricsActor: ActorRef
) extends JobPreparationActor(
      workflowDescriptor = workflowDescriptor,
      jobKey = jobKey,
      factory = null,
      workflowDockerLookupActor = workflowDockerLookupActor,
      initializationData = None,
      serviceRegistryActor = serviceRegistryActor,
      ioActor = ioActor,
      backendSingletonActor = None,
      groupMetricsActor = groupMetricsActor
    ) {

  override private[preparation] lazy val kvStoreKeysToPrefetch = kvStoreKeysForPrefetch

  override private[preparation] lazy val expressionLanguageFunctions = NoIoFunctionSet
  override private[preparation] lazy val dockerHashCredentials = dockerHashCredentialsInput
  override private[preparation] lazy val runtimeAttributeDefinitions = Set.empty
  override private[preparation] lazy val noResponseTimeout = dockerNoResponseTimeoutInput
  override private[preparation] lazy val hasDockerDefinition = true
  override private[preparation] lazy val dockerMirroring = dockerMirroringInput
  override private[preparation] lazy val platform = None

  override private[preparation] def scopedKey(key: String): ScopedKey = scopedKeyMaker.apply(key)
  override private[preparation] def evaluateInputsAndAttributes(valueStore: ValueStore) =
    inputsAndAttributes.getOrElse(super.evaluateInputsAndAttributes(valueStore))

  override private[preparation] def jobExecutionProps(jobDescriptor: BackendJobDescriptor,
                                                      initializationData: Option[BackendInitializationData],
                                                      serviceRegistryActor: ActorRef,
                                                      ioActor: ActorRef,
                                                      backendSingletonActor: Option[ActorRef],
                                                      groupMetricsActor: ActorRef
  ) = Props.empty
}

object JobPreparationTestHelper {
  type ScopedKeyMaker = String => ScopedKey
}
