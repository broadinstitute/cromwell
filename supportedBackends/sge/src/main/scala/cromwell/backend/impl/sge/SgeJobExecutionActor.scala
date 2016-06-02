package cromwell.backend.impl.sge

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend._

import scala.concurrent.{ExecutionContext, Future}

class SgeJobExecutionActor extends BackendJobExecutionActor {
  override protected implicit def ec: ExecutionContext = ???

  override def recover: Future[BackendJobExecutionResponse] = ???

  override def execute: Future[BackendJobExecutionResponse] = ???

  override protected def jobDescriptor: BackendJobDescriptor = ???

  override protected def configurationDescriptor: BackendConfigurationDescriptor = ???

  override def abort: Unit = ???
}
