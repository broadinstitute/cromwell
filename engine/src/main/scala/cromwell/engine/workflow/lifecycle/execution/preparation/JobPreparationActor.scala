package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.{ActorRef, Cancellable, FSM, Props}
import cromwell.backend._
import cromwell.backend.validation.RuntimeAttributesKeys
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.callcaching.docker.DockerHashActor.{DockerHashBackPressure, DockerHashFailureResponse, DockerHashResponseSuccess}
import cromwell.core.callcaching.docker._
import cromwell.core.callcaching._
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation._
import cromwell.engine.workflow.lifecycle.execution.preparation.JobPreparationActor._
import wdl4s._
import wdl4s.values.{WdlString, WdlValue}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class JobPreparationActor(executionData: WorkflowExecutionActorData,
                                     jobKey: BackendJobDescriptorKey,
                                     factory: BackendLifecycleActorFactory,
                                     dockerHashingActor: ActorRef,
                                     initializationData: Option[BackendInitializationData],
                                     serviceRegistryActor: ActorRef,
                                     ioActor: ActorRef,
                                     backendSingletonActor: Option[ActorRef])
  extends FSM[JobPreparationActorState, Option[JobPreparationActorData]] with WorkflowLogging {

  lazy val workflowIdForLogging = workflowDescriptor.id
  
  // Number of attempts to make contact with the docker actor before giving up if no response
  private [preparation] lazy val maxDockerAttempts = 3
  // Amount of time after which the docker request should be considered lost and sent agagin
  private [preparation] lazy val dockerNoResponseTimeout = 30 seconds
  // Amount of time to wait when we get a Backpressure response before sending the request again
  private [preparation] lazy val backpressureWaitTime = 15 seconds
  
  private lazy val workflowDescriptor = executionData.workflowDescriptor
  private val ec = context.system.dispatcher
  private lazy val expressionLanguageFunctions = factory.expressionLanguageFunctions(workflowDescriptor.backendDescriptor, jobKey, initializationData)
  private lazy val dockerHashCredentials = factory.dockerHashCredentials(initializationData)
  
  startWith(Idle, None)

  when(Idle) {
    case Event(Start, _) =>
      evaluateInputsAndAttributes match {
        case Success((inputs, attributes)) => startPreparation(inputs, attributes)
        case Failure(failure) => sendFailureAndStop(failure)
      }
  }

  when(WaitingForDockerHash) {
    case Event(DockerHashResponseSuccess(dockerHash), Some(data)) =>
      handleDockerHashSuccess(dockerHash, data)
    case Event(failureResponse: DockerHashFailureResponse, Some(data)) =>
      log.warning(failureResponse.reason)
      handleDockerHashFailed(failureResponse.reason, data)
    case Event(DockerHashBackPressure(request), Some(data)) =>
      handleDockHashBackPressure(request, data)
    case Event(DockerNoResponseTimeout, Some(data)) =>
      if (data.maxAttemptsReached) handleDockerHashFailed("Docker Hash service is unreachable", data)
      else {
        dockerHashingActor ! data.dockerHashRequest
        stay() using data.newAttempt(newRequestLostTimer)
      }
  }
  
  whenUnhandled {
    case Event(unexpectedMessage, _) =>
      workflowLogger.warn(s"JobPreparation actor received an unexpected message in state $stateName: $unexpectedMessage")
      stay()
  }
  
  private [preparation] def evaluateInputsAndAttributes = {
    for {
      evaluatedInputs <- resolveAndEvaluateInputs(jobKey, workflowDescriptor, expressionLanguageFunctions, executionData.outputStore)
      runtimeAttributes <- prepareRuntimeAttributes(evaluatedInputs)
    } yield (evaluatedInputs, runtimeAttributes)
  }
  
  private def startPreparation(inputs: Map[Declaration, WdlValue], attributes: Map[LocallyQualifiedName, WdlValue]) = {
    def sendDockerRequest(dockerImageId: DockerImageIdentifierWithoutHash) = {
      val dockerHashRequest = DockerHashRequest(dockerImageId, dockerHashCredentials)
      val newData = JobPreparationActorData(dockerHashRequest, inputs, attributes, newRequestLostTimer, maxDockerAttempts)
      dockerHashingActor ! dockerHashRequest
      goto(WaitingForDockerHash) using Option(newData)
    }
    
    def handleDockerValue(value: String) = DockerImageIdentifier.fromString(value) match {
      case Success(dockerImageId: DockerImageIdentifierWithoutHash) => sendDockerRequest(dockerImageId)
      case Success(dockerImageIdWithHash: DockerImageIdentifierWithHash) =>
        // If the docker value already has a hash - we're ok for call caching
        val response = prepareBackendDescriptor(inputs, attributes, CallCachingEligible)
        sendResponseAndStop(response)
      case Failure(failure) => sendFailureAndStop(failure)
    }

    attributes.get(RuntimeAttributesKeys.DockerKey) match {
      case Some(dockerValue) => handleDockerValue(dockerValue.valueString)
      case None =>
        // If there is no docker attribute at all - we're ok for call caching
        val response = prepareBackendDescriptor(inputs, attributes, CallCachingEligible)
        sendResponseAndStop(response)
    }
  }
  
  // Schedule a DockerNoResponseTimeout message to be sent to self
  private final def newRequestLostTimer = context.system.scheduler.scheduleOnce(dockerNoResponseTimeout, self, DockerNoResponseTimeout)(ec, self)
  // Schedule a DockerHashRequest to be sent to the docker hashing actor
  private final def newBackPressureTimer(dockerHashRequest: DockerHashRequest) = context.system.scheduler.scheduleOnce(backpressureWaitTime, dockerHashingActor, dockerHashRequest)(ec, self)

  /*
   * If we get backpressured, recreate both timers
   */
  private def handleDockHashBackPressure(request: DockerHashRequest, data: JobPreparationActorData) = {
    data.cancelTimers()
    val requestLostTimer = newRequestLostTimer
    val backpressureTimer = newBackPressureTimer(data.dockerHashRequest)
    stay() using data.backpressure(requestLostTimer, backpressureTimer)
  }
  
  private def handleDockerHashSuccess(dockerHashResult: DockerHashResult, data: JobPreparationActorData) = {
    data.cancelTimers()
    val dockerValueWithHash = data.dockerHashRequest.dockerImageID.withHash(dockerHashResult).fullName
    val newRuntimeAttributes = data.attributes updated (RuntimeAttributesKeys.DockerKey, WdlString(dockerValueWithHash))
    // We had to ask for the docker hash, which means the attribute had a floating tag, we're NOT ok for call caching
    val response = prepareBackendDescriptor(data.inputs, newRuntimeAttributes, FloatingDockerTagWithHash(dockerValueWithHash))
    sendResponseAndStop(response)
  }

  private def handleDockerHashFailed(failure: String, data: JobPreparationActorData) = {
    data.cancelTimers()
    val response = prepareBackendDescriptor(data.inputs, data.attributes, FloatingDockerTagWithoutHash)
    sendResponseAndStop(response)
  }

  private def sendResponseAndStop(response: CallPreparationActorResponse) = {
    context.parent ! response
    stay()
  }
  
  private def sendFailureAndStop(failure: Throwable) = {
    sendResponseAndStop(CallPreparationFailed(jobKey, failure))
  }
  
  private [preparation] def prepareBackendDescriptor(inputEvaluation: Map[Declaration, WdlValue],
                               runtimeAttributes: Map[LocallyQualifiedName, WdlValue],
                               callCachingEligibility: CallCachingEligibility) = {
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, runtimeAttributes, inputEvaluation, callCachingEligibility)
    BackendJobPreparationSucceeded(jobDescriptor, factory.jobExecutionActorProps(jobDescriptor, initializationData, serviceRegistryActor, ioActor, backendSingletonActor))
  }
  
  private [preparation] def prepareRuntimeAttributes(inputEvaluation: Map[Declaration, WdlValue]): Try[Map[LocallyQualifiedName, WdlValue]] = {
    import RuntimeAttributeDefinition.{addDefaultsToAttributes, evaluateRuntimeAttributes}
    val curriedAddDefaultsToAttributes = addDefaultsToAttributes(factory.runtimeAttributeDefinitions(initializationData), workflowDescriptor.backendDescriptor.workflowOptions) _

    for {
      unevaluatedRuntimeAttributes <- Try(jobKey.call.task.runtimeAttributes)
      evaluatedRuntimeAttributes <- evaluateRuntimeAttributes(unevaluatedRuntimeAttributes, expressionLanguageFunctions, inputEvaluation)
    } yield curriedAddDefaultsToAttributes(evaluatedRuntimeAttributes)
  }
}

object JobPreparationActor {
  case class JobPreparationActorData(dockerHashRequest: DockerHashRequest,
                                     inputs: Map[Declaration, WdlValue],
                                     attributes: Map[LocallyQualifiedName, WdlValue],
                                     requestLostTimer: Cancellable,
                                     maxAttempts: Int,
                                     currentAttempts: Int = 0,
                                     backpressureTimer: Option[Cancellable] = None
                                    ) {
    def cancelTimers(): Unit = {
      requestLostTimer.cancel()
      backpressureTimer foreach { _.cancel() }
      ()
    }
    
    def backpressure(newBackpressureTimer: Cancellable, newRequestLostTimer: Cancellable) = {
      Option(copy(requestLostTimer = newRequestLostTimer, backpressureTimer = Option(newBackpressureTimer)))
    }
    
    def newAttempt(newBackpressureTimer: Cancellable) = {
      Option(copy(requestLostTimer = newBackpressureTimer, currentAttempts = currentAttempts + 1))
    }
    
    def maxAttemptsReached = currentAttempts == maxAttempts
  }

  sealed trait JobPreparationActorState
  case object Idle extends JobPreparationActorState
  case object WaitingForDockerHash extends JobPreparationActorState
  
  private case object DockerNoResponseTimeout
  
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
