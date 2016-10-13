package cromwell.backend.impl.aws

import akka.actor.Props
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}

import scala.concurrent.Future

class AwsJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                           override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  override def execute: Future[BackendJobExecutionResponse] = ???
}

object AwsJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor): Props = Props(new AwsJobExecutionActor(jobDescriptor, configurationDescriptor))
}
