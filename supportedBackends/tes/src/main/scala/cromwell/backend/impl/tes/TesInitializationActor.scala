package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cromwell.backend.standard._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.path.Obsolete._
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilder, PathBuilderFactory}
import wdl4s.TaskCall

import scala.concurrent.Future
import scala.util.Try

case class TesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  tesConfiguration: TesConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = tesConfiguration.configurationDescriptor
}

class TesInitializationActor(params: TesInitializationActorParams)
  extends StandardInitializationActor(params) {

  private val tesConfiguration = params.tesConfiguration

  lazy val pathBuilderFactories: List[PathBuilderFactory] = List(Option(DefaultPathBuilderFactory)).flatten

  override lazy val pathBuilders: List[PathBuilder] =
    pathBuilderFactories map {
      _.withOptions(workflowDescriptor.workflowOptions)(context.system)
    }

  override lazy val workflowPaths: TesWorkflowPaths =
    new TesWorkflowPaths(workflowDescriptor, tesConfiguration.configurationDescriptor.backendConfig, pathBuilders)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    TesRuntimeAttributes.runtimeAttributesBuilder

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.toString)
      File(workflowPaths.workflowRoot).createPermissionedDirectories()
      Option(TesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, tesConfiguration))
    })
  }
}
