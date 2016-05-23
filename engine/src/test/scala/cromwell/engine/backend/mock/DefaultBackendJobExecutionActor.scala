package cromwell.engine.backend.mock

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse}
import cromwell.backend._
import wdl4s.Call
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}

import scala.concurrent.Future

object DefaultBackendJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = Props(DefaultBackendJobExecutionActor(jobDescriptor, configurationDescriptor))
}

case class DefaultBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor, override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {
  override def execute: Future[BackendJobExecutionResponse] = {
    Future.successful(SucceededResponse(jobDescriptor.key, Some(0), (jobDescriptor.call.task.outputs map taskOutputToJobOutput).toMap))
  }
  override def recover = execute

  override def abortJob: Unit = ()
}

class DefaultBackendLifecycleActorFactory(config: Config) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                configurationDescriptor: BackendConfigurationDescriptor): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      configurationDescriptor: BackendConfigurationDescriptor): Props = {
    DefaultBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def workflowFinalizationActorProps(): Option[Props] = None

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           configurationDescriptor: BackendConfigurationDescriptor): WdlStandardLibraryFunctions = NoFunctions
}

