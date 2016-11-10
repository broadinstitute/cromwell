package cromwell.backend.impl.tes

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import cromwell.backend.async.AsyncBackendJobExecutionActor.Execute
import scala.concurrent.{Future, Promise}

final case class TesJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                      override val configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendJobExecutionActor {
  private lazy val completionPromise = Promise[BackendJobExecutionResponse]()
  private var executor: Option[ActorRef] = None

  override def execute: Future[BackendJobExecutionResponse] = {
    for {
      _ <- launchExecutor
      _ = executor foreach { _ ! Execute }
      c <- completionPromise.future
    } yield c
  }

  private def launchExecutor: Future[Unit] = {
    Future {
      val executionProps = TesAsyncBackendJobExecutionActor.props(
        workflowId, jobDescriptor, completionPromise, configurationDescriptor)

      val executorRef = context.actorOf(executionProps, s"TesAsyncBackendJobExecutionActor-$workflowId")
      executor = Option(executorRef)
      ()
    }
  }
}

object TesJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor): Props = {
    Props(new TesJobExecutionActor(jobDescriptor, configurationDescriptor))
  }
}
