package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cromwell.backend.standard._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.path.PathBuilder
import spray.json.JsString
import wom.graph.CommandCallNode

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

case class TesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[CommandCallNode],
  tesConfiguration: TesConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = tesConfiguration.configurationDescriptor
}

class TesInitializationActor(params: TesInitializationActorParams)
  extends StandardInitializationActor(params) {

  private val tesConfiguration = params.tesConfiguration

  override lazy val pathBuilders: Future[List[PathBuilder]] = {
    standardParams.configurationDescriptor.pathBuildersWithDefault(workflowDescriptor.workflowOptions)
  }

  override lazy val workflowPaths: Future[TesWorkflowPaths] = pathBuilders map {
    new TesWorkflowPaths(workflowDescriptor, tesConfiguration.configurationDescriptor.backendConfig, _)
  }

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    TesRuntimeAttributes.runtimeAttributesBuilder(tesConfiguration.runtimeConfig)

  override def validateWorkflowOptions(): Try[Unit] =
    workflowDescriptor.workflowOptions.toMap.get(TesWorkflowOptionKeys.WorkflowExecutionIdentity) match {
      case None => Success(())
      case Some(_: JsString) => Success(())
      case Some(v) => Failure(
        new Exception(s"Workflow option ${TesWorkflowOptionKeys.WorkflowExecutionIdentity} must be a string, was ${v}.")
      )
    }

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    workflowPaths map { paths =>
      publishWorkflowRoot(paths.workflowRoot.toString)
      paths.workflowRoot.createPermissionedDirectories()
      Option(TesBackendInitializationData(paths, runtimeAttributesBuilder, tesConfiguration))
    }
  }
}
