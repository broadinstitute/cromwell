package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.{ActorRef, FSM, Props}
import cromwell.backend._
import cromwell.backend.validation.{DockerValidation, RuntimeAttributesKeys}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.callcaching._
import cromwell.core.logging.WorkflowLogging
import cromwell.docker.DockerHashActor.{DockerHashFailureResponse, DockerHashResponseSuccess}
import cromwell.docker._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation._
import cromwell.engine.workflow.lifecycle.execution.preparation.JobPreparationActor.{DockerNoResponseTimeout, _}
import cromwell.services.keyvalue.KeyValueServiceActor.{KvGet, KvJobKey, KvResponse, ScopedKey}
import wdl4s._
import wdl4s.values.{WdlString, WdlValue}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

class JobPreparationActor(executionData: WorkflowExecutionActorData,
                          jobKey: BackendJobDescriptorKey,
                          factory: BackendLifecycleActorFactory,
                          val dockerHashingActor: ActorRef,
                          initializationData: Option[BackendInitializationData],
                          serviceRegistryActor: ActorRef,
                          ioActor: ActorRef,
                          backendSingletonActor: Option[ActorRef])
  extends FSM[JobPreparationActorState, JobPreparationActorData] with WorkflowLogging with DockerClientHelper {

  override lazy val workflowIdForLogging = workflowDescriptor.id

  // Amount of time after which the docker request should be considered lost and sent agagin
  override protected def backpressureTimeout: FiniteDuration = 10 seconds
  // Amount of time to wait when we get a Backpressure response before sending the request again
  override protected def backpressureRandomizerFactor: Double = 0.5D

  private[preparation] lazy val noResponseTimeout: FiniteDuration = 3 minutes
  
  private lazy val workflowDescriptor = executionData.workflowDescriptor
  private[preparation] lazy val expressionLanguageFunctions = factory.expressionLanguageFunctions(workflowDescriptor.backendDescriptor, jobKey, initializationData)
  private[preparation] lazy val dockerHashCredentials = factory.dockerHashCredentials(initializationData)
  private[preparation] lazy val runtimeAttributeDefinitions = factory.runtimeAttributeDefinitions(initializationData)
  private[preparation] lazy val hasDockerDefinition = runtimeAttributeDefinitions.exists(_.name == DockerValidation.instance.key)

  startWith(Idle, JobPreparationActorNoData)

  context.become(dockerReceive orElse receive)

  when(Idle) {
    case Event(Start, JobPreparationActorNoData) =>
      val keyLookupRequests = kvStoreKeysToPrefetch
      if (keyLookupRequests.nonEmpty) {
        lookupKeyValueEntries(keyLookupRequests)
      } else {
        evaluateInputsAndFetchDockerHashes(KeyValueLookupResults(Map.empty))
      }
  }

  when(FetchingKeyValueStoreEntries) {
    case Event(kvResponse: KvResponse, JobPreparationKeyLookupData(keyLookups)) =>
      keyLookups.withResponse(kvResponse.key, kvResponse) match {
        case newPartialLookup: PartialKeyValueLookups => stay using JobPreparationKeyLookupData(newPartialLookup)
        case finished: KeyValueLookupResults => evaluateInputsAndFetchDockerHashes(finished)
      }
  }

  when(WaitingForDockerHash) {
    case Event(DockerHashResponseSuccess(dockerHash, _), data: JobPreparationHashLookupData) =>
      handleDockerHashSuccess(dockerHash, data)
    case Event(failureResponse: DockerHashFailureResponse, data: JobPreparationHashLookupData) =>
      log.warning(failureResponse.reason)
      handleDockerHashFailed(failureResponse.reason, data)
    case Event(DockerNoResponseTimeout(dockerHashRequest), data: JobPreparationHashLookupData) =>
      handleDockerHashFailed(s"Timed out waiting for a hash for docker image: ${dockerHashRequest.dockerImageID.fullName}", data)
  }

  whenUnhandled {
    case Event(unexpectedMessage, _) =>
      workflowLogger.warn(s"JobPreparation actor received an unexpected message in state $stateName: $unexpectedMessage")
      stay()
  }

  private[preparation] lazy val kvStoreKeysToPrefetch = factory.requestedKeyValueStoreKeys
  private[preparation] def scopedKey(key: String) = ScopedKey(workflowDescriptor.id, KvJobKey(jobKey), key)
  private[preparation] def lookupKeyValueEntries(lookups: Seq[String]) = {
    val keysToLookup: Seq[ScopedKey] = lookups map scopedKey
    keysToLookup foreach { serviceRegistryActor ! KvGet(_) }
    goto(FetchingKeyValueStoreEntries) using JobPreparationKeyLookupData(PartialKeyValueLookups(Map.empty, keysToLookup))
  }

  private [preparation] def evaluateInputsAndAttributes = {
    for {
      evaluatedInputs <- resolveAndEvaluateInputs(jobKey, workflowDescriptor, expressionLanguageFunctions, executionData.outputStore)
      runtimeAttributes <- prepareRuntimeAttributes(evaluatedInputs)
    } yield (evaluatedInputs, runtimeAttributes)
  }

  private def evaluateInputsAndFetchDockerHashes(kvStoreLookupResults: KeyValueLookupResults) = {
    evaluateInputsAndAttributes match {
      case Success((inputs, attributes)) => fetchDockerHashes(kvStoreLookupResults, inputs, attributes)
      case Failure(failure) => sendFailureAndStop(failure)
    }
  }

  private def fetchDockerHashes(kvStoreLookupResults: KeyValueLookupResults, inputs: Map[Declaration, WdlValue], attributes: Map[LocallyQualifiedName, WdlValue]) = {
    def sendDockerRequest(dockerImageId: DockerImageIdentifierWithoutHash) = {
      val dockerHashRequest = DockerHashRequest(dockerImageId, dockerHashCredentials)
      val newData = JobPreparationHashLookupData(kvStoreLookupResults, dockerHashRequest, inputs, attributes)
      sendDockerCommand(dockerHashRequest, noResponseTimeout)
      goto(WaitingForDockerHash) using newData
    }

    def handleDockerValue(value: String) = DockerImageIdentifier.fromString(value) match {
      case Success(dockerImageId: DockerImageIdentifierWithoutHash) if hasDockerDefinition => sendDockerRequest(dockerImageId)
      case Success(_) =>
        // If the docker value already has a hash, or the backend doesn't support docker - 
        // no need to lookup and we're ok for call caching
        val response = prepareBackendDescriptor(inputs, attributes, CallCachingEligible, kvStoreLookupResults.unscoped)
        sendResponseAndStop(response)
      case Failure(failure) => sendFailureAndStop(failure)
    }

    attributes.get(RuntimeAttributesKeys.DockerKey) match {
      case Some(dockerValue) => handleDockerValue(dockerValue.valueString)
      case None =>
        // If there is no docker attribute at all - we're ok for call caching
        val response = prepareBackendDescriptor(inputs, attributes, CallCachingEligible, kvStoreLookupResults.unscoped)
        sendResponseAndStop(response)
    }
  }
  
  private def handleDockerHashSuccess(dockerHashResult: DockerHashResult, data: JobPreparationHashLookupData) = {
    val dockerValueWithHash = data.dockerHashRequest.dockerImageID.withHash(dockerHashResult).fullName
    val newRuntimeAttributes = data.attributes updated (RuntimeAttributesKeys.DockerKey, WdlString(dockerValueWithHash))
    // We had to ask for the docker hash, which means the attribute had a floating tag, we're NOT ok for call caching
    val response = prepareBackendDescriptor(data.inputs, newRuntimeAttributes, FloatingDockerTagWithHash(dockerValueWithHash), data.keyLookupResults.unscoped)
    sendResponseAndStop(response)
  }

  private def handleDockerHashFailed(failure: String, data: JobPreparationHashLookupData) = {
    val response = prepareBackendDescriptor(data.inputs, data.attributes, FloatingDockerTagWithoutHash, data.keyLookupResults.unscoped)
    sendResponseAndStop(response)
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

  private [preparation] def prepareBackendDescriptor(inputEvaluation: Map[Declaration, WdlValue],
                               runtimeAttributes: Map[LocallyQualifiedName, WdlValue],
                               callCachingEligibility: CallCachingEligibility,
                               prefetchedJobStoreEntries: Map[String, KvResponse]): BackendJobPreparationSucceeded = {
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, runtimeAttributes, inputEvaluation, callCachingEligibility, prefetchedJobStoreEntries)
    BackendJobPreparationSucceeded(jobDescriptor, jobExecutionProps(jobDescriptor, initializationData, serviceRegistryActor, ioActor, backendSingletonActor))
  }

  private [preparation] def prepareRuntimeAttributes(inputEvaluation: Map[Declaration, WdlValue]): Try[Map[LocallyQualifiedName, WdlValue]] = {
    import RuntimeAttributeDefinition.{addDefaultsToAttributes, evaluateRuntimeAttributes}
    val curriedAddDefaultsToAttributes = addDefaultsToAttributes(runtimeAttributeDefinitions, workflowDescriptor.backendDescriptor.workflowOptions) _

    for {
      unevaluatedRuntimeAttributes <- Try(jobKey.call.task.runtimeAttributes)
      evaluatedRuntimeAttributes <- evaluateRuntimeAttributes(unevaluatedRuntimeAttributes, expressionLanguageFunctions, inputEvaluation)
    } yield curriedAddDefaultsToAttributes(evaluatedRuntimeAttributes)
  }

  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    message match {
      case request: DockerHashRequest => self ! DockerNoResponseTimeout(request)
      case other => log.warning(s"Unknown request $other timed out")
    }
  }
}

object JobPreparationActor {

  sealed trait JobPreparationActorData
  case object JobPreparationActorNoData extends JobPreparationActorData
  case class JobPreparationKeyLookupData(keyLookups: PartialKeyValueLookups) extends JobPreparationActorData
  private final case class JobPreparationHashLookupData(keyLookupResults: KeyValueLookupResults,
                                                        dockerHashRequest: DockerHashRequest,
                                                        inputs: Map[Declaration, WdlValue],
                                                        attributes: Map[LocallyQualifiedName, WdlValue]) extends JobPreparationActorData

  sealed trait JobPreparationActorState
  case object Idle extends JobPreparationActorState
  case object WaitingForDockerHash extends JobPreparationActorState
  case object FetchingKeyValueStoreEntries extends JobPreparationActorState
  
  private case class DockerNoResponseTimeout(dockerHashRequest: DockerHashRequest)
  
  def props(executionData: WorkflowExecutionActorData,
            jobKey: BackendJobDescriptorKey,
            factory: BackendLifecycleActorFactory,
            dockerHashingActor: ActorRef,
            initializationData: Option[BackendInitializationData],
            serviceRegistryActor: ActorRef,
            ioActor: ActorRef,
            backendSingletonActor: Option[ActorRef]) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new JobPreparationActor(executionData,
      jobKey,
      factory,
      dockerHashingActor,
      initializationData,
      serviceRegistryActor,
      ioActor,
      backendSingletonActor)).withDispatcher(EngineDispatcher)
  }
}
