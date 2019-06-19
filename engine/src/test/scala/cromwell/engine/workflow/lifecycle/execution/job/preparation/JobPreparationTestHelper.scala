package cromwell.engine.workflow.lifecycle.execution.job.preparation

import wdl.draft2.model.LocallyQualifiedName
import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.TestProbe
import cromwell.backend._
import cromwell.core.WorkflowId
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData
import cromwell.engine.workflow.lifecycle.execution.job.preparation.JobPreparationTestHelper._
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.services.keyvalue.KeyValueServiceActor.{KvJobKey, ScopedKey}
import common.validation.ErrorOr.ErrorOr
import org.specs2.mock.Mockito
import wom.expression.NoIoFunctionSet
import wom.graph.{CommandCallNode, WomIdentifier}
import wom.values.{WomEvaluatedCallInputs, WomValue}

import scala.concurrent.duration.FiniteDuration

class JobPreparationTestHelper(implicit val system: ActorSystem) extends Mockito {
  val executionData = mock[WorkflowExecutionActorData]
  val workflowDescriptor = mock[EngineWorkflowDescriptor]
  val backendDescriptor = mock[BackendWorkflowDescriptor]
  val workflowId = WorkflowId.randomId()
  workflowDescriptor.backendDescriptor returns backendDescriptor
  workflowDescriptor.id returns workflowId
  workflowDescriptor.possiblyNotRootWorkflowId returns workflowId.toPossiblyNotRoot
  workflowDescriptor.rootWorkflowId returns workflowId.toRoot
  workflowDescriptor.rootWorkflow returns workflowDescriptor
  executionData.workflowDescriptor returns workflowDescriptor
  val jobKey = mock[BackendJobDescriptorKey]
  val call = CommandCallNode(WomIdentifier("JobPreparationSpec_call"), null, null, null, Set.empty, null, None)
  jobKey.call returns call
  jobKey.node returns call
  jobKey.index returns None
  jobKey.attempt returns 1
  val serviceRegistryProbe = TestProbe()
  val ioActor = TestProbe()
  val workflowDockerLookupActor = TestProbe()

  val scopedKeyMaker: ScopedKeyMaker = key => ScopedKey(workflowId, KvJobKey("correct.horse.battery.staple", None, 1), key)

  def buildTestJobPreparationActor(backpressureTimeout: FiniteDuration,
                                   noResponseTimeout: FiniteDuration,
                                   dockerHashCredentials: List[Any],
                                   inputsAndAttributes: ErrorOr[(WomEvaluatedCallInputs, Map[LocallyQualifiedName, WomValue])],
                                   kvStoreKeysForPrefetch: List[String]) = {

    Props(new TestJobPreparationActor(
      kvStoreKeysForPrefetch = kvStoreKeysForPrefetch,
      dockerHashCredentialsInput = dockerHashCredentials,
      backpressureWaitTimeInput = backpressureTimeout,
      dockerNoResponseTimeoutInput = noResponseTimeout,
      inputsAndAttributes = inputsAndAttributes,
      workflowDescriptor = workflowDescriptor,
      jobKey = jobKey,
      workflowDockerLookupActor = workflowDockerLookupActor.ref,
      serviceRegistryActor = serviceRegistryProbe.ref,
      ioActor = ioActor.ref,
      scopedKeyMaker))
  }
}

private[preparation] class TestJobPreparationActor(kvStoreKeysForPrefetch: List[String],
                                                   dockerHashCredentialsInput: List[Any],
                                                   backpressureWaitTimeInput: FiniteDuration,
                                                   dockerNoResponseTimeoutInput: FiniteDuration,
                                                   inputsAndAttributes: ErrorOr[(WomEvaluatedCallInputs, Map[LocallyQualifiedName, WomValue])],
                                                   workflowDescriptor: EngineWorkflowDescriptor,
                                                   jobKey: BackendJobDescriptorKey,
                                                   workflowDockerLookupActor: ActorRef,
                                                   serviceRegistryActor: ActorRef,
                                                   ioActor: ActorRef,
                                                   scopedKeyMaker: ScopedKeyMaker) extends JobPreparationActor(workflowDescriptor = workflowDescriptor,
                                                                                                  jobKey = jobKey,
                                                                                                  factory = null,
                                                                                                  workflowDockerLookupActor = workflowDockerLookupActor,
                                                                                                  initializationData = None,
                                                                                                  serviceRegistryActor = serviceRegistryActor,
                                                                                                  ioActor = ioActor,
                                                                                                  backendSingletonActor = None) {

  override lazy val kvStoreKeysToPrefetch = kvStoreKeysForPrefetch

  override private[preparation] lazy val expressionLanguageFunctions = NoIoFunctionSet
  override private[preparation] lazy val dockerHashCredentials = dockerHashCredentialsInput
  override private[preparation] lazy val noResponseTimeout = dockerNoResponseTimeoutInput
  override private[preparation] lazy val hasDockerDefinition = true

  override def scopedKey(key: String) = scopedKeyMaker.apply(key)
  override def evaluateInputsAndAttributes(valueStore: ValueStore) = inputsAndAttributes

  override private[preparation] def jobExecutionProps(jobDescriptor: BackendJobDescriptor,
                                                      initializationData: Option[BackendInitializationData],
                                                      serviceRegistryActor: ActorRef,
                                                      ioActor: ActorRef,
                                                      backendSingletonActor: Option[ActorRef]) = Props.empty
}

object JobPreparationTestHelper {
  type ScopedKeyMaker = String => ScopedKey
}
