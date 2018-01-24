package cromwell.backend.impl.bcs

import akka.actor.ActorRef
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.oss.OssPathBuilderFactory
import net.ceedubs.ficus.Ficus._
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._

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
  private implicit val system = context.system

  lazy val ossPathBuilderFactory: Option[OssPathBuilderFactory] = {
    for {
      endpoint <- configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.endpoint")
      accessId <- configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.access-id")
      accessKey <- configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.access-key")
      securityToken = configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.security-token")
    } yield new OssPathBuilderFactory(endpoint, accessId, accessKey, securityToken)
  }


  override lazy val pathBuilders: Future[List[PathBuilder]] =
    ossPathBuilderFactory.toList.traverse(_.withOptions(workflowDescriptor.workflowOptions)).map(_ ++ Option(DefaultPathBuilder))


  override lazy val workflowPaths: Future[BcsWorkflowPaths] = pathBuilders map {
    BcsWorkflowPaths(workflowDescriptor, bcsConfiguration.configurationDescriptor.backendConfig, _)
  }

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    BcsRuntimeAttributes.runtimeAttributesBuilder(bcsConfiguration.runtimeConfig)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    pathBuilders map { builders => BcsMount.pathBuilders = builders}

    for {
      paths <- workflowPaths
      builders <- pathBuilders
    } yield Option(BcsBackendInitializationData(paths, runtimeAttributesBuilder, bcsConfiguration, builders))
  }
}
