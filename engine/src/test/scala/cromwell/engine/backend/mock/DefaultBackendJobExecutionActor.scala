package cromwell.engine.backend.mock

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobSucceededResponse}
import cromwell.backend._
import cromwell.core.{CallOutputs, NoIoFunctionSet}
import wom.expression.IoFunctionSet
import wom.graph.TaskCallNode

import scala.concurrent.{ExecutionContext, Future}

object DefaultBackendJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor) = Props(DefaultBackendJobExecutionActor(jobDescriptor, configurationDescriptor))
}

case class DefaultBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor, override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {
  override def execute: Future[BackendJobExecutionResponse] = {
    Future.successful(JobSucceededResponse(jobDescriptor.key, Some(0), CallOutputs((jobDescriptor.taskCall.outputPorts map taskOutputToJobOutput).toMap), None, Seq.empty, dockerImageUsed = None))
  }

  override def recover = execute

  override def abort(): Unit = ()
}

class DefaultBackendLifecycleActorFactory(val name: String, val configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                ioActor: ActorRef,
                                                calls: Set[TaskCallNode],
                                                serviceRegistryActor: ActorRef,
                                                restarting: Boolean): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      ioActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    DefaultBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData],
                                           ioActorEndpoint: ActorRef,
                                           ec: ExecutionContext): IoFunctionSet = NoIoFunctionSet
}

