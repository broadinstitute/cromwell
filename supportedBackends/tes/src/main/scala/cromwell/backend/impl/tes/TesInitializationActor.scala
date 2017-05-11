package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.standard._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilder, PathBuilderFactory}
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import wdl4s.TaskCall
import net.ceedubs.ficus.Ficus._

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

  /**
    * If the backend sets a gcs authentication mode, try to create a PathBuilderFactory with it.
    */
  lazy val gcsPathBuilderFactory: Option[GcsPathBuilderFactory] = {
    configurationDescriptor.backendConfig.as[Option[String]]("filesystems.gcs.auth") map { configAuth =>
      val googleConfiguration = GoogleConfiguration(configurationDescriptor.globalConfig)
      googleConfiguration.auth(configAuth) match {
        case Valid(auth) => GcsPathBuilderFactory(auth, googleConfiguration.applicationName)
        case Invalid(error) => throw new MessageAggregation {
          override def exceptionContext: String = "Failed to parse gcs auth configuration"

          override def errorMessages: Traversable[String] = error.toList
        }
      }
    }
  }

  lazy val pathBuilderFactories: List[PathBuilderFactory] =
    List(gcsPathBuilderFactory, Option(DefaultPathBuilderFactory)).flatten

  override lazy val pathBuilders: List[PathBuilder] =
    pathBuilderFactories map { _.withOptions(workflowDescriptor.workflowOptions)(context.system) }

  override lazy val workflowPaths: TesWorkflowPaths =
    new TesWorkflowPaths(workflowDescriptor, tesConfiguration.configurationDescriptor.backendConfig, pathBuilders)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    TesRuntimeAttributes.runtimeAttributesBuilder(tesConfiguration.runtimeConfig)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.toString)
      workflowPaths.workflowRoot.createPermissionedDirectories()
      Option(TesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, tesConfiguration))
    })
  }
}
