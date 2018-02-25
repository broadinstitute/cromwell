package cromwell.engine.backend.mock

import akka.actor.Props
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobFailedNonRetryableResponse, JobFailedRetryableResponse}

import scala.concurrent.Future

object RetryableBackendJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = Props(RetryableBackendJobExecutionActor(jobDescriptor, configurationDescriptor))
}

final case class RetryableBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor, override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  val attempts = 3

  override def execute: Future[BackendJobExecutionResponse] = {
    if (jobDescriptor.key.attempt < attempts) {
      Future.successful(JobFailedRetryableResponse(jobDescriptor.key, new RuntimeException("An apparent transient Exception!"), None))
    }
    else {
      Future.successful(JobFailedNonRetryableResponse(jobDescriptor.key, new RuntimeException("A permanent Exception! Yikes, what a pickle!"), None))
    }
  }

  override def recover = execute

  override def abort(): Unit = ()
}
