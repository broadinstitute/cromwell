package cromwell.engine.workflow

import akka.actor.{ActorSystem, Props}
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestKitSpec
import cromwell.core.{HogGroup, SimpleIoActor, WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse, WorkflowDescriptorMaterializationResult}

import scala.concurrent.Await

trait WorkflowDescriptorBuilderForSpecs {

  implicit val awaitTimeout = CromwellTestKitSpec.TimeoutDuration
  implicit val actorSystem: ActorSystem
  lazy val ioActor = actorSystem.actorOf(SimpleIoActor.props)

  def createMaterializedEngineWorkflowDescriptor(id: WorkflowId, workflowSources: WorkflowSourceFilesCollection): EngineWorkflowDescriptor = {
    import akka.pattern.ask
    implicit val timeout = akka.util.Timeout(awaitTimeout)
    implicit val ec = actorSystem.dispatcher

    val serviceRegistryIgnorer = actorSystem.actorOf(Props.empty)
    val actor = actorSystem.actorOf(MaterializeWorkflowDescriptorActor.props(serviceRegistryIgnorer, id, importLocalFilesystem = false, ioActorProxy = ioActor, hogGroup = HogGroup("testcase")), "MaterializeWorkflowDescriptorActor-" + id.id)
    val workflowDescriptorFuture = actor.ask(
      MaterializeWorkflowDescriptorCommand(workflowSources, ConfigFactory.load)
    ).mapTo[WorkflowDescriptorMaterializationResult]

    Await.result(workflowDescriptorFuture map {
      case MaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor) => workflowDescriptor
      case MaterializeWorkflowDescriptorFailureResponse(reason) => throw reason
    }, awaitTimeout)
  }
}
