package cromwell.engine.backend.mock

import akka.actor.Props
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionFailedRetryableResponse, BackendJobExecutionResponse, BackendJobExecutionSucceededResponse}

import scala.concurrent.Future

object RetryableBackendJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = Props(RetryableBackendJobExecutionActor(jobDescriptor, configurationDescriptor))
}

case class RetryableBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor, override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  val attempts = 3

  override def execute: Future[BackendJobExecutionResponse] = {
    if(jobDescriptor.key.attempt < attempts)
      Future.successful(BackendJobExecutionFailedRetryableResponse(jobDescriptor.key, new RuntimeException("An apparent transient Exception!")))
    else
      Future.successful(BackendJobExecutionSucceededResponse(jobDescriptor.key, (jobDescriptor.call.task.outputs map taskOutputToJobOutput).toMap))
  }
  override def recover = execute

  override def abortJob: Unit = ()
}
