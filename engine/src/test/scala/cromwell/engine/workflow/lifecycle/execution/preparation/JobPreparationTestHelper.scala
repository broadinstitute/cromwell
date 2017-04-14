package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.TestProbe
import cromwell.backend._
import cromwell.core.WorkflowId
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData
import cromwell.services.keyvalue.KeyValueServiceActor.{KvJobKey, ScopedKey}
import org.specs2.mock.Mockito
import wdl4s._
import wdl4s.values.WdlValue

import scala.concurrent.duration.FiniteDuration
import scala.util.Try
import JobPreparationTestHelper._
import wdl4s.expression.NoFunctions

class JobPreparationTestHelper(implicit val system: ActorSystem) extends Mockito {
  val executionData = mock[WorkflowExecutionActorData]
  val workflowDescriptor = mock[EngineWorkflowDescriptor]
  val backendDescriptor = mock[BackendWorkflowDescriptor]
  workflowDescriptor.backendDescriptor returns backendDescriptor
  executionData.workflowDescriptor returns workflowDescriptor
  val jobKey = mock[BackendJobDescriptorKey]
  val serviceRegistryProbe = TestProbe()
  val ioActor = TestProbe()
  val dockerHashingActor = TestProbe()

  val workflowId = WorkflowId.randomId()
  val scopedKeyMaker: ScopedKeyMaker = key => ScopedKey(workflowId, KvJobKey("correct.horse.battery.staple", None, 1), key)

  def buildTestJobPreparationActor(backpressureTimeout: FiniteDuration,
                                   noResponseTimeout: FiniteDuration,
                                   dockerHashCredentials: List[Any],
                                   inputsAndAttributes: Try[(Map[Declaration, WdlValue], Map[wdl4s.LocallyQualifiedName, WdlValue])],
                                   kvStoreKeysForPrefetch: List[String]) = {

    Props(new TestJobPreparationActor(
      kvStoreKeysForPrefetch = kvStoreKeysForPrefetch,
      dockerHashCredentialsInput = dockerHashCredentials,
      backpressureWaitTimeInput = backpressureTimeout,
      dockerNoResponseTimeoutInput = noResponseTimeout,
      inputsAndAttributes = inputsAndAttributes,
      executionData = executionData,
      jobKey = jobKey,
      dockerHashingActor = dockerHashingActor.ref,
      serviceRegistryActor = serviceRegistryProbe.ref,
      ioActor = ioActor.ref,
      scopedKeyMaker))
  }
}

private[preparation] class TestJobPreparationActor(kvStoreKeysForPrefetch: List[String],
                                                   dockerHashCredentialsInput: List[Any],
                                                   backpressureWaitTimeInput: FiniteDuration,
                                                   dockerNoResponseTimeoutInput: FiniteDuration,
                                                   inputsAndAttributes: Try[(Map[Declaration, WdlValue], Map[wdl4s.LocallyQualifiedName, WdlValue])],
                                                   executionData: WorkflowExecutionActorData,
                                                   jobKey: BackendJobDescriptorKey,
                                                   dockerHashingActor: ActorRef,
                                                   serviceRegistryActor: ActorRef,
                                                   ioActor: ActorRef,
                                                   scopedKeyMaker: ScopedKeyMaker) extends JobPreparationActor(executionData = executionData,
                                                                                                  jobKey = jobKey,
                                                                                                  factory = null,
                                                                                                  dockerHashingActor = dockerHashingActor,
                                                                                                  initializationData = None,
                                                                                                  serviceRegistryActor = serviceRegistryActor,
                                                                                                  ioActor = ioActor,
                                                                                                  backendSingletonActor = None) {

  override lazy val kvStoreKeysToPrefetch = kvStoreKeysForPrefetch

  override private[preparation] lazy val expressionLanguageFunctions = NoFunctions
  override private[preparation] lazy val dockerHashCredentials = dockerHashCredentialsInput
  override private[preparation] lazy val noResponseTimeout = dockerNoResponseTimeoutInput
  override private[preparation] lazy val hasDockerDefinition = true
  override protected def backpressureTimeout: FiniteDuration = backpressureWaitTimeInput

  override def scopedKey(key: String) = scopedKeyMaker.apply(key)
  override def evaluateInputsAndAttributes = inputsAndAttributes

  override private[preparation] def jobExecutionProps(jobDescriptor: BackendJobDescriptor,
                                                      initializationData: Option[BackendInitializationData],
                                                      serviceRegistryActor: ActorRef,
                                                      ioActor: ActorRef,
                                                      backendSingletonActor: Option[ActorRef]) = Props.empty
}

object JobPreparationTestHelper {
  type ScopedKeyMaker = String => ScopedKey
}
