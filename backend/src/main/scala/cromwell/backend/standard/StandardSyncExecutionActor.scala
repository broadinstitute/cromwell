package cromwell.backend.standard

import akka.actor.SupervisorStrategy.{Decider, Stop}
import akka.actor.{ActorRef, OneForOneStrategy, Props}
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, Recover}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import cromwell.core.Dispatcher
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{Future, Promise}

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

  context.become(startup orElse super.receive)

  private def startup: Receive = {
    case AbortJobCommand =>
      context.parent ! AbortedResponse(jobDescriptor.key)
      context.stop(self)
  }

  private def running(executor: ActorRef): Receive = {
    case AbortJobCommand =>
      executor ! AbortJobCommand
    case abortResponse: AbortedResponse =>
      context.parent ! abortResponse
      context.stop(self)
    case KvPair(key, Some(jobId)) if key.key == jobIdKey =>
      // Successful operation ID lookup during recover.
      executor ! Recover(StandardAsyncJob(jobId))
    case KvKeyLookupFailed(_) =>
      // Missed operation ID lookup during recover, fall back to execute.
      executor ! Execute
    case KvFailure(_, e) =>
      // Failed operation ID lookup during recover, crash and let the supervisor deal with it.
      completionPromise.tryFailure(e)
      throw new RuntimeException(s"Failure attempting to look up job id for key ${jobDescriptor.key}", e)
  }

  /**
    * This "synchronous" actor isn't finished until this promise finishes over in the asynchronous version.
    *
    * Still not sure why the AsyncBackendJobExecutionActor doesn't wait for an Akka message instead of using Scala promises.
    */
  lazy val completionPromise: Promise[BackendJobExecutionResponse] = Promise[BackendJobExecutionResponse]()

  override def execute: Future[BackendJobExecutionResponse] = {
    val executorRef = createAsyncRef()
    context.become(running(executorRef) orElse super.receive)
    executorRef ! Execute
    completionPromise.future
  }

  override def recover: Future[BackendJobExecutionResponse] = {
    val executorRef = createAsyncRef()
    context.become(running(executorRef) orElse super.receive)
    val kvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt)
    val kvGet = KvGet(ScopedKey(jobDescriptor.workflowDescriptor.id, kvJobKey, jobIdKey))
    serviceRegistryActor ! kvGet
    completionPromise.future
  }

  def createAsyncParams(): StandardAsyncExecutionActorParams = {
    DefaultStandardAsyncExecutionActorParams(
      standardParams.jobIdKey,
      standardParams.serviceRegistryActor,
      standardParams.jobDescriptor,
      standardParams.configurationDescriptor,
      standardParams.backendInitializationDataOption,
      completionPromise
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
    throw new NotImplementedError("Abort is implemented via a custom receive of the message AbortJobCommand.")
  }

  // Supervision strategy: if the async actor throws an exception, stop the actor and fail the job.
  def jobFailingDecider: Decider = {
    case exception: Exception =>
      completionPromise.tryFailure(
        new RuntimeException(s"${createAsyncRefName()} failed and didn't catch its exception.", exception))
      Stop
  }

  override val supervisorStrategy: OneForOneStrategy = OneForOneStrategy()(jobFailingDecider)
}
