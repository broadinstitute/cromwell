package cromwell.engine.workflow

import akka.actor.{ActorSystem, Props}
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse, WorkflowDescriptorMaterializationResult}
import cromwell.engine.{EngineWorkflowDescriptor, WorkflowSourceFiles}

import scala.concurrent.Await

trait WorkflowDescriptorBuilder {

  implicit val awaitTimeout = CromwellTestkitSpec.timeoutDuration
  implicit val actorSystem: ActorSystem

  def createMaterializedEngineWorkflowDescriptor(id: WorkflowId, workflowSources: WorkflowSourceFiles): EngineWorkflowDescriptor = {
    import akka.pattern.ask
    implicit val timeout = akka.util.Timeout(awaitTimeout)
    implicit val ec = actorSystem.dispatcher

    val serviceRegistryIgnorer = actorSystem.actorOf(Props.empty)
    val actor = actorSystem.actorOf(MaterializeWorkflowDescriptorActor.props(serviceRegistryIgnorer, id), "MaterializeWorkflowDescriptorActor-" + id.id)
    val workflowDescriptorFuture = actor.ask(
      MaterializeWorkflowDescriptorCommand(workflowSources, ConfigFactory.load)
    ).mapTo[WorkflowDescriptorMaterializationResult]

    Await.result(workflowDescriptorFuture map {
      case MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor) => workflowDescriptor
      case MaterializeWorkflowDescriptorFailureResponse(reason) => throw reason
    }, awaitTimeout)
  }
}
