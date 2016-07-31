package cromwell.backend.sfs

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.Execute
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}

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
  * @param jobDescriptor The job to execute.
  * @param configurationDescriptor The configuration.
  * @param asyncPropsCreator A function that can create the specific asynchronous backend.
  */
class SharedFileSystemJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                        override val configurationDescriptor: BackendConfigurationDescriptor,
                                        asyncPropsCreator: Promise[BackendJobExecutionResponse] => Props)
  extends BackendJobExecutionActor {

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
  }

  /**
    * This "synchronous" actor isn't finished until this promise finishes over in the asynchronous version.
    *
    * Still not sure why the AsyncBackendJobExecutionActor doesn't wait for an Akka message instead of using Scala promises.
    */
  private val completionPromise = Promise[BackendJobExecutionResponse]()

  override def execute: Future[BackendJobExecutionResponse] = {
    val executorRef = context.actorOf(asyncPropsCreator(completionPromise), "SharedFileSystemAsyncJobExecutionActor")
    context.become(running(executorRef) orElse super.receive)
    executorRef ! Execute
    completionPromise.future
  }

  override def abort() = {
    throw new NotImplementedError("Abort is implemented via a custom receive of the message AbortJobCommand.")
  }
}
