package cromwell.engine

import java.nio.file.Paths

import akka.actor.FSM.NullFunction
import akka.actor._
import akka.event.Logging
import com.google.api.client.util.ExponentialBackOff
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendActor.Prepare
import cromwell.backend.{BackendActor, DefaultBackendFactory}
import cromwell.backend.config.BackendConfiguration
import cromwell.backend.model._
import cromwell.engine
import wdl4s._
import wdl4s.values.{WdlString, WdlValue}
import cromwell.engine.CallActor.{CallActorData, CallActorState}
import cromwell.engine.backend._
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.slick.Execution
import cromwell.engine.workflow.WorkflowActor.{InitialStartCall, UseCachedCall}
import cromwell.engine.workflow.{CallKey, WorkflowActor}
import cromwell.instrumentation.Instrumentation.Monitor
import cromwell.logging.WorkflowLogger

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

import scala.util.{Try, Success, Failure}


object CallActor {

  sealed trait CallActorMessage
  case object Initialize extends CallActorMessage
  sealed trait StartMode extends CallActorMessage
  case object Start extends StartMode
//  final case class Resume(jobKey: JobKey) extends StartMode { override val executionMessage = CallExecutionActor.Resume(jobKey) }
  final case class UseCachedCall(cachedBackendCall: BackendCall, backendCall: BackendCall) extends StartMode
  final case class RegisterCallAbortFunction(abortFunction: AbortFunction) extends CallActorMessage
  case object AbortCall extends CallActorMessage
  final case class ExecutionFinished(call: Call, executionResult: cromwell.engine.backend.ExecutionResult) extends CallActorMessage

  sealed trait CallActorState
  case object CallNotStarted extends CallActorState
  case object CallRunningAbortUnavailable extends CallActorState
  case object CallRunningAbortAvailable extends CallActorState
  case object CallRunningAbortRequested extends CallActorState
  case object CallAborting extends CallActorState
  case object CallDone extends CallActorState

  /** The `WorkflowActor` will drop `TerminalCallMessage`s that it is unable to immediately process.  This message
    * represents the acknowledgment from the `WorkflowActor` that the `callMessage` was processed and the status
    * change has been successfully committed to the database. */
  final case class Ack(callMessage: WorkflowActor.TerminalCallMessage) extends CallActorMessage

  /** Message to self to retry sending the `callMessage` to the `WorkflowActor` in the absence of an `Ack`. */
  final case class Retry(callMessage: WorkflowActor.CallMessage) extends CallActorMessage

  /** FSM data class with an optional abort function and exponential backoff. */
  case class CallActorData(abortFunction: Option[AbortFunction] = None,
                           backoff: Option[ExponentialBackOff] = None,
                           timer: Option[Cancellable] = None) {

    /**
      * Schedule a single retry of sending the specified message to the specified actor.  Return a copy of this
      * `CallActorData` holding a reference to the timer to allow for its cancellation in the event an `Ack` is
      * received prior to its expiration.
      */
    def copyWithRetry(actor: Actor, callMessage: WorkflowActor.CallMessage)(implicit ec: ExecutionContext): CallActorData = {
      val timer = actor.context.system.scheduler.scheduleOnce(backoff.get.nextBackOffMillis().millis) {
        actor.self ! Retry(callMessage)
      }
      this.copy(timer = Option(timer))
    }

    /** Attempt to cancel the last timer if there is one registered. */
    def cancelTimer() = timer foreach { _.cancel() }
  }

  val CallCounter = Monitor.minMaxCounter("calls-running")

  def props(key: CallKey, locallyQualifiedInputs: CallInputs, workflowDescriptor: WorkflowDescriptor): Props =
    Props(new CallActor(key, locallyQualifiedInputs, workflowDescriptor))
}

/** Actor to manage the execution of a single call. */
class CallActor(key: CallKey, locallyQualifiedInputs: CallInputs, workflowDescriptor: WorkflowDescriptor)
  extends LoggingFSM[CallActorState, CallActorData] with CromwellActor {

  import CallActor._

  type CallOutputs = Map[String, WdlValue]

  startWith(CallNotStarted, CallActorData())
  CallCounter.increment()

  val actorSystem = context.system
  implicit val ec = actorSystem.dispatcher

  // TODO: [gaurav] Need this to not be a var and move it to somewhere else
  var backend: Backend = _

  val call = key.scope
  val akkaLogger = Logging(context.system, classOf[CallActor])
  val logger = WorkflowLogger(
    "CallActor",
    workflowDescriptor,
    akkaLogger = Option(akkaLogger),
    callTag = Option(key.tag)
  )

  // Called on every state transition.
  onTransition {
    case _ -> CallDone => shutDown()
    case fromState -> toState =>
      // Only log this at debug - these states are never seen or used outside of the CallActor itself.
      logger.debug(s"transitioning from $fromState to $toState.")
  }

  when(CallNotStarted) {
    case Event(_@Initialize, _) =>
      sendStartMessage()
      stay()
    case Event(startMode: StartMode, _) =>
      // There's no special Retry/Ack handling required for CallStarted message, the WorkflowActor can always
      // handle those immediately.
      context.parent ! WorkflowActor.CallStarted(key)
//      val backendCall = backend.bindCall(workflowDescriptor, key, locallyQualifiedInputs, AbortRegistrationFunction(registerAbortFunction))
      //TODO: Take care of the 'Backend' type
      val backendActor: BackendActor = createBackendActor()
      subscribe(self, backendActor)
//      val executionActorName = s"CallExecutionActor-${workflowDescriptor.id}-${call.unqualifiedName}"
//      context.actorOf(CallExecutionActor.props(), executionActorName) ! backendActor
      goto(CallRunningAbortUnavailable)
    case Event(AbortCall, _) => handleFinished(call, AbortedExecution)
  }

  when(CallRunningAbortUnavailable) {
    case Event(RegisterCallAbortFunction(newAbortFunction), data) =>
      val updatedData = data.copy(abortFunction = Option(newAbortFunction))
      goto(CallRunningAbortAvailable) using updatedData
    case Event(AbortCall, _) => goto(CallRunningAbortRequested)
  }

  when(CallRunningAbortAvailable) {
    case Event(RegisterCallAbortFunction(newAbortFunction), data) =>
      logger.warn("An existing abort function was overwritten with a new one. This might indicate unexpected state transitions in the backend.")
      val updatedData = data.copy(abortFunction = Option(newAbortFunction))
      goto(CallRunningAbortAvailable) using updatedData
    case Event(AbortCall, data) => tryAbort(data)
  }

  when(CallRunningAbortRequested) {
    case Event(RegisterCallAbortFunction(abortFunction), data) =>
      val updatedData = data.copy(abortFunction = Option(abortFunction))
      tryAbort(updatedData)
  }

  // When in CallAborting, the only message being listened for is the CallComplete message (which is already
  // handled by whenUnhandled.)
  when(CallAborting) { NullFunction }

  when(CallDone) {
    case Event(e, _) =>
      logger.error(s"received unexpected event $e while in state $stateName")
      stay()
  }

  import cromwell.backend.model._

  whenUnhandled {
    case Event(Subscribed, _) =>
      sender ! Prepare
      stay()
    case Event(TaskStatus(Status.Created, _), _) =>
      sender ! cromwell.backend.BackendActor.Execute
      stay()
    case Event(TaskStatus(Status.Canceled | Status.Failed, throwable:Throwable), _) =>
      val failureHandle = new FailedExecution(throwable)
      self ! ExecutionFinished(call, failureHandle)
      stay
    case Event(TaskStatus(Status.Succeeded, result: SuccesfulResult), _) =>
      val callOps: Map[String, CallOutput] = result.outputs.map(output => output._1 -> CallOutput(output._2, None))
      val successResult = SuccessfulExecution(callOps, Seq.empty, 0, new ExecutionHash(this.hashCode().toString, None))
      self ! ExecutionFinished(call, successResult)
      stay()

    case Event(ExecutionFinished(finishedCall, executionResult), _) => handleFinished(finishedCall, executionResult)
    case Event(Retry(message), data) =>
      logger.debug(s"CallActor retrying ${message.callKey.tag}")
      // Continue retrying this message until it is Acked.
      context.parent ! message
      val updatedData = data.copyWithRetry(this, message)
      stay() using updatedData
    case Event(Ack(message: WorkflowActor.TerminalCallMessage), data) =>
      logger.debug(s"CallActor received ack for ${message.callKey.tag}")
      data.cancelTimer()
      val updatedData = data.copy(backoff = None, timer = None)
      goto(CallDone) using updatedData
    case Event(e, _) =>
      logger.error(s"received unhandled event $e while in state $stateName")
      stay()
  }

  private def tryAbort(data: CallActorData): CallActor.this.State = {
    data.abortFunction match {
      case Some(af) =>
        logger.info("Abort function called.")
        af.function()
        goto(CallAborting) using data
      case None =>
        logger.warn("Call abort failed because the provided abort function was null.")
        goto(CallRunningAbortRequested) using data
    }
  }

  private def handleFinished(call: Call, executionResult: cromwell.engine.backend.ExecutionResult): State = {

    def createBackoff: Option[ExponentialBackOff] = Option(
      new ExponentialBackOff.Builder()
        .setInitialIntervalMillis(100)
        .setMaxIntervalMillis(3000)
        .setMaxElapsedTimeMillis(Integer.MAX_VALUE)
        .setMultiplier(1.1)
        .build()
    )

    val message = executionResult match {
      case SuccessfulExecution(outputs, executionEvents, returnCode, hash, resultsClonedFrom) =>
        WorkflowActor.CallCompleted(key, outputs, executionEvents, returnCode, if (workflowDescriptor.writeToCache) Option(hash) else None, resultsClonedFrom)
      case AbortedExecution => WorkflowActor.CallAborted(key)
      case FailedExecution(e, returnCode) =>
        logger.error("Failing call: " + e.getMessage, e)
        WorkflowActor.CallFailed(key, returnCode, e.getMessage)
    }

    context.parent ! message
    // The WorkflowActor will drop TerminalCallMessages that it is unable to immediately process.
    // Retry sending this message to the WorkflowActor until it is Acked.
    val updatedData = stateData.copy(backoff = createBackoff).copyWithRetry(this, message)
    stay using updatedData
  }

  private def shutDown(): Unit = {
    logger.debug(s"done, shutting down.")
    CallCounter.decrement()
    context.stop(self)
  }

  private def registerAbortFunction(abortFunction: AbortFunction): Unit = {
    self ! CallActor.RegisterCallAbortFunction(abortFunction)
  }

  //----------------------------------------------------------------------------------------//
  //                                        NEW CHANGES                                     //
  //----------------------------------------------------------------------------------------//

  private def subscribe(subscriberActorRef: ActorRef, backendActor: BackendActor) = {
    val subscription = Subscription[ActorRef](new ExecutionEvent, subscriberActorRef)
    backendActor.self ! BackendActor.SubscribeToEvent[ActorRef](subscription)
  }

  // Currently this method will create a backend
  def createBackendActor(): BackendActor = {
    val backendConfig = BackendConfiguration.apply()
    val backend = backendConfig.getDefaultBackend()
    val task = buildTaskDescriptor()
    val backendActor: BackendActor = DefaultBackendFactory.getBackend(backend.initClass, actorSystem, task)
    backendActor
  }

  /**
    * Returns a TaskDescriptor, given the context of a task
    * @return
    */
  private def buildTaskDescriptor(): TaskDescriptor = {
    val name = call.fullyQualifiedName
    val user = System.getProperty("user.name")

    val engineFunctions = new DefaultWorkflowEngineFunctions(workflowDescriptor.ioManager, workflowDescriptor.wfContext)
    val commandTry = call.instantiateCommandLine(locallyQualifiedInputs, engineFunctions)
    val command = commandTry match {
      case Success(command) => command
      case Failure(error) => throw new IllegalStateException(s"Failed to instantiate command! Message = ${error.getMessage}", error)
    }
    val commandArgs: Seq[String] = command.split(" ")
    val runtimeAttributes = call.task.runtimeAttributes.attrs.map{case (k,v) => (k,v.head)}
    TaskDescriptor(name, user, command, Option(commandArgs), s"${workflowDescriptor.name}-${workflowDescriptor.id}" , locallyQualifiedInputs, call.task.outputs, runtimeAttributes)
  }

  private def wdlStringToScalaString(input: WdlValue): String = {
    // 'toWdlString' returns a string wrapped in double quotes. We need to strip that off
    input.toWdlString.substring(1, input.toWdlString.length - 1)
  }

  private def sendStartMessage() = {
    def registerAbortFunction(abortFunction: AbortFunction): Unit = {}
    val backendCall = backend.bindCall(workflowDescriptor, key, locallyQualifiedInputs, AbortRegistrationFunction(registerAbortFunction))
    val log = backendCall.workflowLoggerWithCall("WorkflowActor", Option(akkaLogger))

    def loadCachedBackendCallAndMessage(descriptor: WorkflowDescriptor, cachedExecution: Execution) = {
      descriptor.namespace.resolve(cachedExecution.callFqn) match {
        case Some(c: Call) =>
          val cachedCall = backend.bindCall(
            descriptor,
            CallKey(c, Option(cachedExecution.index)),
            locallyQualifiedInputs,
            AbortRegistrationFunction(registerAbortFunction)
          )
          log.info(s"Call Caching: Cache hit. Using UUID(${cachedCall.workflowDescriptor.shortId}):${cachedCall.key.tag} as results for UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.key.tag}")
          context.parent ! WorkflowActor.UseCachedCall(key, CallActor.UseCachedCall(cachedCall, backendCall))
        case _ =>
          log.error(s"Call Caching: error when resolving '${cachedExecution.callFqn}' in workflow with execution ID ${cachedExecution.workflowExecutionId}: falling back to normal execution")
          context.parent ! InitialStartCall(key, CallActor.Start)
      }
    }

    /* Tries to use the cached Execution to send a UseCachedCall message.  If anything fails, send an InitialStartCall message */
    def loadCachedCallOrInitiateCall(cachedDescriptor: Try[WorkflowDescriptor], cachedExecution: Execution) = cachedDescriptor match {
      case Success(descriptor) => loadCachedBackendCallAndMessage(descriptor, cachedExecution)
      case Failure(ex) =>
        log.error(s"Call Caching: error when loading workflow with execution ID ${cachedExecution.workflowExecutionId}: falling back to normal execution", ex)
        context.parent ! InitialStartCall(key, CallActor.Start)
    }

    if (backendCall.workflowDescriptor.readFromCache) {
      backendCall.hash map { hash =>
        globalDataAccess.getExecutionsWithResuableResultsByHash(hash.overallHash) onComplete {
          case Success(executions) if executions.nonEmpty =>
            val cachedExecution = executions.head
            globalDataAccess.getWorkflow(cachedExecution.workflowExecutionId) onComplete { cachedDescriptor =>
              loadCachedCallOrInitiateCall(cachedDescriptor, cachedExecution)
            }
          case Success(_) =>
            log.info(s"Call Caching: cache miss")
            context.parent ! InitialStartCall(key, CallActor.Start)
          case Failure(ex) =>
            log.error(s"Call Caching: Failed to look up executions that matched hash '$hash'. Falling back to normal execution", ex)
            context.parent ! InitialStartCall(key, CallActor.Start)
        }
      } recover { case e =>
        log.error(s"Failed to calculate hash for call '${backendCall.key.tag}'.", e)
        val failedExecution = FailedExecution(e)
        handleFinished(call, failedExecution)
        //        context.parent ! WorkflowActor.CallFailed(callKey)
        //        scheduleTransition(WorkflowFailed)
      }
    }
    else {
      log.info(s"Call caching 'readFromCache' is turned off, starting call")
      context.parent ! InitialStartCall(key, CallActor.Start)
    }
  }
}
