package cromwell.engine.workflow.lifecycle.execution.job.preparation

import _root_.wdl.draft2.model._
import akka.actor.{ActorRef, FSM, Props}
import cats.data.Validated.{Invalid, Valid}
import common.exception.MessageAggregation
import common.validation.ErrorOr.ErrorOr
import cromwell.backend._
import cromwell.backend.validation.DockerValidation
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.callcaching._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.{Dispatcher, DockerConfiguration}
import cromwell.docker.DockerInfoActor.{DockerInfoSuccessResponse, DockerInformation, DockerSize}
import cromwell.docker._
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.WorkflowDockerLookupActor.{WorkflowDockerLookupFailure, WorkflowDockerTerminalFailure}
import cromwell.engine.workflow.lifecycle.execution.CallMetadataHelper
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation._
import cromwell.engine.workflow.lifecycle.execution.job.preparation.JobPreparationActor._
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{CallMetadataKeys, MetadataEvent, MetadataValue}
import wom.RuntimeAttributesKeys
import wom.callable.Callable.InputDefinition
import wom.values._

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

/**
  * Prepares a job for the backend. The goal of this actor is to build a BackendJobDescriptor.
  * Those are the steps:
  * 1) Evaluate call inputs and runtime attributes.
  * 3) If a docker hash needs to be looked up, do so. Otherwise go to step 5.
  * 4) When the result comes back, infer the right callCachingEligibility.
  * 5) If Key/Value pairs need to be looked up, do so. Otherwise go to step 6.
  * 6) Build the JobDescriptor, respond and stop.
  */
class JobPreparationActor(workflowDescriptor: EngineWorkflowDescriptor,
                          jobKey: BackendJobDescriptorKey,
                          factory: BackendLifecycleActorFactory,
                          val workflowDockerLookupActor: ActorRef,
                          initializationData: Option[BackendInitializationData],
                          val serviceRegistryActor: ActorRef,
                          ioActor: ActorRef,
                          backendSingletonActor: Option[ActorRef])
  extends FSM[JobPreparationActorState, JobPreparationActorData] with WorkflowLogging with CallMetadataHelper {

  override lazy val workflowIdForLogging = workflowDescriptor.possiblyNotRootWorkflowId
  override lazy val workflowIdForCallMetadata = workflowDescriptor.possiblyNotRootWorkflowId
  override lazy val rootWorkflowIdForLogging = workflowDescriptor.rootWorkflowId

  private[preparation] lazy val noResponseTimeout: FiniteDuration = 3 minutes
  private[preparation] val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  private[preparation] lazy val expressionLanguageFunctions = factory.expressionLanguageFunctions(workflowDescriptor.backendDescriptor, jobKey, initializationData, ioActor, ioEc)
  private[preparation] lazy val dockerHashCredentials = factory.dockerHashCredentials(workflowDescriptor.backendDescriptor, initializationData)
  private[preparation] lazy val runtimeAttributeDefinitions = factory.runtimeAttributeDefinitions(initializationData)
  private[preparation] lazy val hasDockerDefinition = runtimeAttributeDefinitions.exists(_.name == DockerValidation.instance.key)

  startWith(Idle, JobPreparationActorNoData)

  when(Idle) {
    case Event(Start(valueStore), JobPreparationActorNoData) =>
      evaluateInputsAndAttributes(valueStore) match {
        case Valid((inputs, attributes)) => fetchDockerHashesIfNecessary(inputs, attributes)
        case Invalid(failure) => sendFailureAndStop(new MessageAggregation {
          override def exceptionContext: String = s"Call input and runtime attributes evaluation failed for ${jobKey.call.localName}"
          override def errorMessages: Traversable[String] = failure.toList
        })
      }
  }

  when(WaitingForDockerHash) {
    case Event(DockerInfoSuccessResponse(DockerInformation(dockerHash, dockerSize), _), data: JobPreparationDockerLookupData) =>
      handleDockerHashSuccess(dockerHash, dockerSize, data)
    case Event(WorkflowDockerLookupFailure(reason, _), data: JobPreparationDockerLookupData) =>
      workflowLogger.warn("Docker lookup failed", reason)
      handleDockerHashFailed(data)
    case Event(WorkflowDockerTerminalFailure(reason, _), _: JobPreparationDockerLookupData) =>
      sendFailureAndStop(reason)
  }

  when(FetchingKeyValueStoreEntries) {
    case Event(kvResponse: KvResponse, data @ JobPreparationKeyLookupData(keyLookups, maybeCallCachingEligible, dockerSize, inputs, attributes)) =>
      keyLookups.withResponse(kvResponse.key, kvResponse) match {
        case newPartialLookup: PartialKeyValueLookups => stay using data.copy(keyLookups = newPartialLookup)
        case finished: KeyValueLookupResults => 
          sendResponseAndStop(prepareBackendDescriptor(inputs, attributes, maybeCallCachingEligible, finished.unscoped, dockerSize))
      }
  }

  whenUnhandled {
    case Event(unexpectedMessage, _) =>
      workflowLogger.warn(s"JobPreparation actor received an unexpected message in state $stateName: $unexpectedMessage")
      stay()
  }

  private[preparation] lazy val kvStoreKeysToPrefetch: Seq[String] = factory.requestedKeyValueStoreKeys ++ factory.defaultKeyValueStoreKeys
  private[preparation] def scopedKey(key: String) = ScopedKey(workflowDescriptor.id, KvJobKey(jobKey), key)
  private[preparation] def lookupKeyValueEntries(inputs: WomEvaluatedCallInputs,
                                                 attributes: Map[LocallyQualifiedName, WomValue],
                                                 maybeCallCachingEligible: MaybeCallCachingEligible,
                                                 dockerSize: Option[DockerSize]) = {
    val keysToLookup = kvStoreKeysToPrefetch map scopedKey
    keysToLookup foreach { serviceRegistryActor ! KvGet(_) }
    goto(FetchingKeyValueStoreEntries) using JobPreparationKeyLookupData(PartialKeyValueLookups(Map.empty, keysToLookup), maybeCallCachingEligible, dockerSize, inputs, attributes)
  }

  private [preparation] def evaluateInputsAndAttributes(valueStore: ValueStore): ErrorOr[(WomEvaluatedCallInputs, Map[LocallyQualifiedName, WomValue])] = {
    import common.validation.ErrorOr.ShortCircuitingFlatMap
    for {
      evaluatedInputs <- resolveAndEvaluateInputs(jobKey, expressionLanguageFunctions, valueStore)
      runtimeAttributes <- prepareRuntimeAttributes(evaluatedInputs)
    } yield (evaluatedInputs, runtimeAttributes)
  }

  private def fetchDockerHashesIfNecessary(inputs: WomEvaluatedCallInputs, attributes: Map[LocallyQualifiedName, WomValue]) = {
    def sendDockerRequest(dockerImageId: DockerImageIdentifier) = {
      val dockerHashRequest = DockerInfoRequest(dockerImageId, dockerHashCredentials)
      val newData = JobPreparationDockerLookupData(dockerHashRequest, inputs, attributes)
      workflowDockerLookupActor ! dockerHashRequest
      goto(WaitingForDockerHash) using newData
    }

    def handleDockerValue(value: String) = DockerImageIdentifier.fromString(value) match {
        // If the backend supports docker, lookup is enabled, and we got a tag - we need to lookup the hash
      case Success(dockerImageId: DockerImageIdentifierWithoutHash) if hasDockerDefinition && DockerConfiguration.instance.enabled => 
        sendDockerRequest(dockerImageId)
        // If the backend supports docker, we got a tag but lookup is disabled, continue with no call caching and no hash
      case Success(dockerImageId: DockerImageIdentifierWithoutHash) if hasDockerDefinition =>
        lookupKvsOrBuildDescriptorAndStop(inputs, attributes, FloatingDockerTagWithoutHash(dockerImageId.fullName), None)

      // If the backend doesn't support docker - no need to lookup and we're ok for call caching
      case Success(_: DockerImageIdentifierWithoutHash) if !hasDockerDefinition => 
        lookupKvsOrBuildDescriptorAndStop(inputs, attributes, NoDocker, None)

      // Even if the docker image has a hash, we need to (try to) find out the size, so send a request
      case Success(dockerImageId: DockerImageIdentifierWithHash) =>
        sendDockerRequest(dockerImageId)

      case Failure(failure) => sendFailureAndStop(failure)
    }

    attributes.get(RuntimeAttributesKeys.DockerKey) match {
      case Some(dockerValue) => handleDockerValue(dockerValue.valueString)
      case None =>
        // If there is no docker attribute at all - we're ok for call caching
        lookupKvsOrBuildDescriptorAndStop(inputs, attributes, NoDocker, None)
    }
  }

  /**
    * Either look up KVs or build the JobDescriptor with empty KV map and respond.
    */
  private def lookupKvsOrBuildDescriptorAndStop(inputs: WomEvaluatedCallInputs,
                                                attributes: Map[LocallyQualifiedName, WomValue],
                                                maybeCallCachingEligible: MaybeCallCachingEligible,
                                                dockerSize: Option[DockerSize]) = {
    if (kvStoreKeysToPrefetch.nonEmpty) lookupKeyValueEntries(inputs, attributes, maybeCallCachingEligible, dockerSize)
    else sendResponseAndStop(prepareBackendDescriptor(inputs, attributes, maybeCallCachingEligible, Map.empty, dockerSize))
  }

  private def handleDockerHashSuccess(dockerHashResult: DockerHashResult, dockerSize: Option[DockerSize], data: JobPreparationDockerLookupData) = {
    val hashValue = data.dockerHashRequest.dockerImageID match {
      case withoutHash: DockerImageIdentifierWithoutHash => withoutHash.withHash(dockerHashResult)
      case withHash => withHash
    }
    dockerSize.foreach(sendCompressedDockerSizeToMetadata)
    lookupKvsOrBuildDescriptorAndStop(data.inputs, data.attributes, DockerWithHash(hashValue.fullName), dockerSize)
  }
  
  private def sendCompressedDockerSizeToMetadata(dockerSize: DockerSize) = {
    val event = MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.CompressedDockerSize), MetadataValue(dockerSize.compressedSize))
    serviceRegistryActor ! PutMetadataAction(event)
  }

  private def handleDockerHashFailed(data: JobPreparationDockerLookupData) = {
    val floatingDockerTag = data.dockerHashRequest.dockerImageID.fullName
    lookupKvsOrBuildDescriptorAndStop(data.inputs, data.attributes, FloatingDockerTagWithoutHash(floatingDockerTag), None)
  }

  private def sendResponseAndStop(response: CallPreparationActorResponse) = {
    context.parent ! response
    stay()
  }

  private def sendFailureAndStop(failure: Throwable) = {
    sendResponseAndStop(CallPreparationFailed(jobKey, failure))
  }

  // 'jobExecutionProps' is broken into a separate function for TestJobPreparationActor to override:
  private[preparation] def jobExecutionProps(jobDescriptor: BackendJobDescriptor,
                                             initializationData: Option[BackendInitializationData],
                                             serviceRegistryActor: ActorRef,
                                             ioActor: ActorRef,
                                             backendSingletonActor: Option[ActorRef]) = factory.jobExecutionActorProps(jobDescriptor, initializationData, serviceRegistryActor, ioActor, backendSingletonActor)

  private[preparation] def prepareBackendDescriptor(inputEvaluation: WomEvaluatedCallInputs,
                                                    runtimeAttributes: Map[LocallyQualifiedName, WomValue],
                                                    maybeCallCachingEligible: MaybeCallCachingEligible,
                                                    prefetchedJobStoreEntries: Map[String, KvResponse],
                                                    dockerSize: Option[DockerSize]): BackendJobPreparationSucceeded = {
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, runtimeAttributes, inputEvaluation, maybeCallCachingEligible, dockerSize, prefetchedJobStoreEntries)
    BackendJobPreparationSucceeded(jobDescriptor, jobExecutionProps(jobDescriptor, initializationData, serviceRegistryActor, ioActor, backendSingletonActor))
  }

  private [preparation] def prepareRuntimeAttributes(inputEvaluation: Map[InputDefinition, WomValue]): ErrorOr[Map[LocallyQualifiedName, WomValue]] = {
    import RuntimeAttributeDefinition.{addDefaultsToAttributes, evaluateRuntimeAttributes}
    val curriedAddDefaultsToAttributes = addDefaultsToAttributes(runtimeAttributeDefinitions, workflowDescriptor.backendDescriptor.workflowOptions) _
    val unevaluatedRuntimeAttributes = jobKey.call.callable.runtimeAttributes
    evaluateRuntimeAttributes(unevaluatedRuntimeAttributes, expressionLanguageFunctions, inputEvaluation) map curriedAddDefaultsToAttributes
  }
}

object JobPreparationActor {

  sealed trait JobPreparationActorData
  case object JobPreparationActorNoData extends JobPreparationActorData
  private final case class JobPreparationKeyLookupData(keyLookups: PartialKeyValueLookups,
                                                       maybeCallCachingEligible: MaybeCallCachingEligible,
                                                       dockerSize: Option[DockerSize],
                                                       inputs: WomEvaluatedCallInputs,
                                                       attributes: Map[LocallyQualifiedName, WomValue]) extends JobPreparationActorData
  private final case class JobPreparationDockerLookupData(dockerHashRequest: DockerInfoRequest,
                                                          inputs: WomEvaluatedCallInputs,
                                                          attributes: Map[LocallyQualifiedName, WomValue]) extends JobPreparationActorData

  sealed trait JobPreparationActorState
  case object Idle extends JobPreparationActorState
  case object WaitingForDockerHash extends JobPreparationActorState
  case object FetchingKeyValueStoreEntries extends JobPreparationActorState

  def props(workflowDescriptor: EngineWorkflowDescriptor,
            jobKey: BackendJobDescriptorKey,
            factory: BackendLifecycleActorFactory,
            workflowDockerLookupActor: ActorRef,
            initializationData: Option[BackendInitializationData],
            serviceRegistryActor: ActorRef,
            ioActor: ActorRef,
            backendSingletonActor: Option[ActorRef]) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new JobPreparationActor(workflowDescriptor,
      jobKey,
      factory,
      workflowDockerLookupActor = workflowDockerLookupActor,
      initializationData,
      serviceRegistryActor = serviceRegistryActor,
      ioActor = ioActor,
      backendSingletonActor = backendSingletonActor)).withDispatcher(EngineDispatcher)
  }
}
