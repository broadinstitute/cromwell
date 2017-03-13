package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cromwell.backend.standard._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
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

  override lazy val workflowPaths: TesWorkflowPaths =
    new TesWorkflowPaths(workflowDescriptor, tesConfiguration.configurationDescriptor.backendConfig, pathBuilders)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    TesRuntimeAttributes.runtimeAttributesBuilder(tesConfiguration.runtimeAttrsConfig)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.toString)
      workflowPaths.workflowRoot.createPermissionedDirectories()
      Option(TesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, tesConfiguration))
    })
  }
}
