package cromwell.backend.sfs

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.Execute
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}

import scala.concurrent.{Future, Promise}

class SharedFileSystemJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                        override val configurationDescriptor: BackendConfigurationDescriptor,
                                        asyncPropsCreator: Promise[BackendJobExecutionResponse] => Props)
  extends BackendJobExecutionActor {

    context.become(sharedStartup orElse super.receive)

    def sharedStartup: Receive = {
      case AbortJobCommand =>
        context.parent ! AbortedResponse(jobDescriptor.key)
        context.stop(self)
      case abortResponse: AbortedResponse =>
        context.parent ! abortResponse
        context.stop(self)
    }

    def sharedRunning(executor: ActorRef): Receive = {
      case AbortJobCommand =>
        executor ! AbortJobCommand
      case abortResponse: AbortedResponse =>
        context.parent ! abortResponse
        context.stop(self)
    }

    override def execute: Future[BackendJobExecutionResponse] = {
      // Still not sure why the AsyncJobExecutionActor doesn't wait for an Akka message instead of using Scala promises
      val completionPromise = Promise[BackendJobExecutionResponse]()
      val executorRef = context.actorOf(asyncPropsCreator(completionPromise), "SharedFileSystemAsyncJobExecutionActor")
      context.become(sharedRunning(executorRef) orElse super.receive)
      executorRef ! Execute
      completionPromise.future
    }
  }
