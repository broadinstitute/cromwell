package cromwell.backend.standard

import akka.actor.SupervisorStrategy.{Decider, Stop}
import akka.actor.{ActorRef, OneForOneStrategy, Props}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobAbortedResponse, JobNotFoundException}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, Reconnect, ReconnectToAbort, Recover}
import cromwell.core.Dispatcher
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{Future, Promise}
import scala.util.control.NoStackTrace

trait StandardSyncExecutionActorParams extends StandardJobExecutionActorParams {
  /** The class for creating an async backend. */
  def asyncJobExecutionActorClass: Class[_ <: StandardAsyncExecutionActor]
}

case class DefaultStandardSyncExecutionActorParams
(
  override val jobIdKey: String,
  override val serviceRegistryActor: ActorRef,
  override val ioActor: ActorRef,
  override val jobDescriptor: BackendJobDescriptor,
  override val configurationDescriptor: BackendConfigurationDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val backendSingletonActorOption: Option[ActorRef],
  override val asyncJobExecutionActorClass: Class[_ <: StandardAsyncExecutionActor],
  override val minimumRuntimeSettings: MinimumRuntimeSettings
) extends StandardSyncExecutionActorParams

/**
  * Facade to the asynchronous execution actor.
  *
  * Creates the asynchronous execution actor, then relays messages to that actor.
  *
  * NOTE: Although some methods return futures due to the (current) contract in BJEA/ABJEA, this actor only executes
  * during the receive, and does not launch new runnables/futures from inside "receive".
  *
  * Thus there are no vars, and the context switches during "receive", once the asynchronous actor has been created.
  *
  * The backend synchronous actor is a proxy for an asynchronous actor that does the actual heavy lifting.
  *
  * - Synchronous params contain the class of the asynchronous actor.
  * - Synchronous actor creates `Promise` that it will wait for.
  * - Synchronous actor passes the `Promise` into a `StandardAsyncExecutionActorParams`.
  * - Synchronous actor creates a `Props` using the asynchronous class plus asynchronous params.
  * - Synchronous actor creates an asynchronous actor using the `Props`.
  * - Synchronous actor waits for the `Promise` to complete.
  * - Asynchronous actor runs.
  * - Asynchronous actor completes the promise with a success or failure.
  */
class StandardSyncExecutionActor(val standardParams: StandardSyncExecutionActorParams)
  extends BackendJobExecutionActor {

  override val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor
  override val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor
  val jobIdKey: String = standardParams.jobIdKey
  val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor

  context.become(startup orElse receive)

  private def startup: Receive = {
    case AbortJobCommand =>
      context.parent ! JobAbortedResponse(jobDescriptor.key)
      context.stop(self)
  }

  private def running(executor: ActorRef): Receive = {
    case AbortJobCommand =>
      executor ! AbortJobCommand
    case abortResponse: JobAbortedResponse =>
      context.parent ! abortResponse
      context.stop(self)
    case KvFailure(_, e) =>
      // Failed operation ID lookup, crash and let the supervisor deal with it.
      completionPromise.tryFailure(e)
      throw new RuntimeException(s"Failure attempting to look up job id for key ${jobDescriptor.key}", e)
  }
  
  private def recovering(executor: ActorRef): Receive = running(executor).orElse {
    case KvPair(key, jobId) if key.key == jobIdKey =>
      // Successful operation ID lookup.
      executor ! Recover(StandardAsyncJob(jobId))
    case KvKeyLookupFailed(_) =>
      // Missed operation ID lookup, fall back to execute.
      executor ! Execute
  }

  private def reconnectingToAbort(executor: ActorRef): Receive = running(executor).orElse {
    case KvPair(key, jobId) if key.key == jobIdKey =>
      // Successful operation ID lookup.
      executor ! ReconnectToAbort(StandardAsyncJob(jobId))
    case KvKeyLookupFailed(_) =>
      // Can't find an operation ID for this job, respond with a JobNotFound and don't start the job
      completionPromise.tryFailure(JobNotFoundException(jobDescriptor.key))
      ()
  }

  private def reconnecting(executor: ActorRef): Receive = running(executor).orElse {
    case KvPair(key, jobId) if key.key == jobIdKey =>
      // Successful operation ID lookup.
      executor ! Reconnect(StandardAsyncJob(jobId))
    case KvKeyLookupFailed(_) =>
      // Can't find an operation ID for this job, respond with a JobNotFound and don't start the job
      completionPromise.tryFailure(JobNotFoundException(jobDescriptor.key))
      ()
  }

  /**
    * This "synchronous" actor isn't finished until this promise finishes over in the asynchronous version.
    *
    * Still not sure why the AsyncBackendJobExecutionActor doesn't wait for an Akka message instead of using Scala promises.
    */
  lazy val completionPromise: Promise[BackendJobExecutionResponse] = Promise[BackendJobExecutionResponse]()

  override def execute: Future[BackendJobExecutionResponse] = {
    val executorRef = createAsyncRef()
    context.become(running(executorRef) orElse receive)
    executorRef ! Execute
    completionPromise.future
  }

  private def onRestart(nextReceive: ActorRef => Receive) = {
    val executorRef = createAsyncRef()
    context.become(nextReceive(executorRef) orElse receive)
    val kvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt)
    val kvGet = KvGet(ScopedKey(jobDescriptor.workflowDescriptor.id, kvJobKey, jobIdKey))
    serviceRegistryActor ! kvGet
    completionPromise.future
  }
  
  override def recover: Future[BackendJobExecutionResponse] = {
    onRestart(recovering)
  }
  
  override def reconnectToAborting: Future[BackendJobExecutionResponse] = {
    onRestart(reconnectingToAbort)
  }

  override def reconnect: Future[BackendJobExecutionResponse] = {
    onRestart(reconnecting)
  }

  def createAsyncParams(): StandardAsyncExecutionActorParams = {
    DefaultStandardAsyncExecutionActorParams(
      standardParams.jobIdKey,
      standardParams.serviceRegistryActor,
      standardParams.ioActor,
      standardParams.jobDescriptor,
      standardParams.configurationDescriptor,
      standardParams.backendInitializationDataOption,
      standardParams.backendSingletonActorOption,
      completionPromise,
      standardParams.minimumRuntimeSettings
    )
  }

  def createAsyncProps(): Props = {
    val asyncParams = createAsyncParams()
    Props(standardParams.asyncJobExecutionActorClass, asyncParams)
  }

  def createAsyncRefName(): String = {
    standardParams.asyncJobExecutionActorClass.getSimpleName
  }

  def createAsyncRef(): ActorRef = {
    val props = createAsyncProps().withDispatcher(Dispatcher.BackendDispatcher)
    val name = createAsyncRefName()
    context.actorOf(props, name)
  }

  override def abort(): Unit = {
    throw new UnsupportedOperationException("Abort is implemented via a custom receive of the message AbortJobCommand.")
  }

  // Supervision strategy: if the async actor throws an exception, stop the actor and fail the job.
  def jobFailingDecider: Decider = {
    case exception: Exception =>
      completionPromise.tryFailure(
        new RuntimeException(s"${createAsyncRefName()} failed and didn't catch its exception. This condition has been handled and the job will be marked as failed.", exception) with NoStackTrace)
      Stop
  }

  override val supervisorStrategy: OneForOneStrategy = OneForOneStrategy()(jobFailingDecider)
}
