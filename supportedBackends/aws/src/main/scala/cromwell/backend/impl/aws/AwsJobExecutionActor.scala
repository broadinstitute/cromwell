package cromwell.backend.impl.aws

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, Recover}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import cromwell.core.Dispatcher
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{Future, Promise}

class AwsJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                           override val configurationDescriptor: BackendConfigurationDescriptor,
                           serviceRegistryActor: ActorRef) extends BackendJobExecutionActor {

  // Copypasta

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
    case KvPair(key, id@Some(jobId)) if key.key == jobIdKey =>
      // Successful operation ID lookup during recover.
      executor ! recoverMessage(jobId)
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
    */
  private lazy val completionPromise = Promise[BackendJobExecutionResponse]()

  override def execute: Future[BackendJobExecutionResponse] = {
    val executorRef = context.actorOf(asyncProps, "SharedFileSystemAsyncJobExecutionActor")
    context.become(running(executorRef) orElse super.receive)
    executorRef ! Execute
    completionPromise.future
  }

  override def recover: Future[BackendJobExecutionResponse] = {
    val executorRef = context.actorOf(asyncProps, "SharedFileSystemAsyncJobExecutionActor")
    context.become(running(executorRef) orElse super.receive)
    val kvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt)
    val kvGet = KvGet(ScopedKey(jobDescriptor.workflowDescriptor.id, kvJobKey, jobIdKey))
    serviceRegistryActor ! kvGet
    completionPromise.future
  }

  override def abort() = {
    throw new NotImplementedError("Abort is implemented via a custom receive of the message AbortJobCommand.")
  }

  // aws specific

  private val jobIdKey = AwsJobId.JobIdKey

  private def recoverMessage(jobId: String): Recover = {
    Recover(AwsJobId(jobId))
  }

  private lazy val asyncProps: Props = {
    Props(
      new AwsAsyncJobExecutionActor(jobDescriptor, completionPromise, configurationDescriptor, serviceRegistryActor)
    ).withDispatcher(Dispatcher.BackendDispatcher)
  }
}

object AwsJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor,
            serviceRegistryActor: ActorRef): Props = {
    Props(
      new AwsJobExecutionActor(jobDescriptor, configurationDescriptor, serviceRegistryActor)
    ).withDispatcher(Dispatcher.BackendDispatcher)
  }
}
