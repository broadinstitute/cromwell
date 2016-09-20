package cromwell.engine.backend.mock

import akka.actor.Props
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import cromwell.backend.BackendJobExecutionActor.{FailedRetryableResponse, BackendJobExecutionResponse, SucceededResponse}

import scala.concurrent.Future

object RetryableBackendJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = Props(RetryableBackendJobExecutionActor(jobDescriptor, configurationDescriptor))
}

case class RetryableBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor, override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  val attempts = 3

  override def execute: Future[BackendJobExecutionResponse] = {
    if (jobDescriptor.key.attempt < attempts)
      Future.successful(FailedRetryableResponse(jobDescriptor.key, new RuntimeException("An apparent transient Exception!"), None))
    else
      Future.successful(SucceededResponse(jobDescriptor.key, Some(0), (jobDescriptor.call.task.outputs map taskOutputToJobOutput).toMap, None, Seq.empty))
  }

  override def recover = execute

  override def abort(): Unit = ()
}
