package cromwell.backend.impl.bcs

import akka.actor.ActorRef
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.path.PathBuilder

import scala.concurrent.Future
import wom.graph.CommandCallNode


final case class BcsInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[CommandCallNode],
  bcsConfiguration: BcsConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = bcsConfiguration.configurationDescriptor
}

final class BcsInitializationActor(params: BcsInitializationActorParams)
  extends StandardInitializationActor(params) {

  private val bcsConfiguration = params.bcsConfiguration

  override lazy val pathBuilders: Future[List[PathBuilder]] =
    standardParams.configurationDescriptor.pathBuildersWithDefault(workflowDescriptor.workflowOptions)

  override lazy val workflowPaths: Future[BcsWorkflowPaths] = pathBuilders map {
    BcsWorkflowPaths(workflowDescriptor, bcsConfiguration.configurationDescriptor.backendConfig, _)
  }

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder = BcsRuntimeAttributes.runtimeAttributesBuilder(bcsConfiguration.runtimeConfig)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    pathBuilders map { builders => BcsMount.pathBuilders = builders}

    for {
      paths <- workflowPaths
      builders <- pathBuilders
    } yield Option(BcsBackendInitializationData(paths, runtimeAttributesBuilder, bcsConfiguration, builders))

  }

}
