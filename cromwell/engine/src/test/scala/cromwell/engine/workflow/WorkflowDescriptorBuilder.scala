package cromwell.engine.workflow

import akka.actor.{ActorSystem, Props}
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestKitSpec
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse, WorkflowDescriptorMaterializationResult}

import scala.concurrent.Await

trait WorkflowDescriptorBuilder {

  implicit val awaitTimeout = CromwellTestKitSpec.TimeoutDuration
  implicit val actorSystem: ActorSystem

  def createMaterializedEngineWorkflowDescriptor(id: WorkflowId, workflowSources: WorkflowSourceFilesCollection): EngineWorkflowDescriptor = {
    import akka.pattern.ask
    implicit val timeout = akka.util.Timeout(awaitTimeout)
    implicit val ec = actorSystem.dispatcher

    val serviceRegistryIgnorer = actorSystem.actorOf(Props.empty)
    val actor = actorSystem.actorOf(MaterializeWorkflowDescriptorActor.props(serviceRegistryIgnorer, id, importLocalFilesystem = false), "MaterializeWorkflowDescriptorActor-" + id.id)
    val workflowDescriptorFuture = actor.ask(
      MaterializeWorkflowDescriptorCommand(workflowSources, ConfigFactory.load)
    ).mapTo[WorkflowDescriptorMaterializationResult]

    Await.result(workflowDescriptorFuture map {
      case MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor) => workflowDescriptor
      case MaterializeWorkflowDescriptorFailureResponse(reason) => throw reason
    }, awaitTimeout)
  }
}
