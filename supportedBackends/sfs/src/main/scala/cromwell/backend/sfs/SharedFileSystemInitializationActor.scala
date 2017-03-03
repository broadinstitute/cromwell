package cromwell.backend.sfs

import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.BackendInitializationData
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.standard.{StandardExpressionFunctions, StandardInitializationActor, StandardInitializationActorParams}
import cromwell.backend.wfs.WorkflowPathBuilder
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilder, PathBuilderFactory}
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._

import scala.concurrent.Future
import scala.util.Try

class SharedFileSystemInitializationActor(standardParams: StandardInitializationActorParams)
  extends StandardInitializationActor(standardParams) {

  /**
    * If the backend sets a gcs authentication mode, try to create a PathBuilderFactory with it.
    */
  lazy val gcsPathBuilderFactory: Option[GcsPathBuilderFactory] = {
    configurationDescriptor.backendConfig.as[Option[String]]("filesystems.gcs.auth") map { configAuth =>
      GoogleConfiguration(configurationDescriptor.globalConfig).auth(configAuth) match {
        case Valid(auth) => GcsPathBuilderFactory(auth)
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

  override lazy val workflowPaths: WorkflowPaths =
    WorkflowPathBuilder.workflowPaths(configurationDescriptor, workflowDescriptor, pathBuilders)

  override lazy val expressionFunctions: Class[_ <: StandardExpressionFunctions] =
    classOf[SharedFileSystemExpressionFunctions]

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.pathAsString)
      workflowPaths.workflowRoot.createPermissionedDirectories()
      Option(initializationData)
    })
  }
}
