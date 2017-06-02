package cromwell.backend.sfs

import cats.data.Validated.{Invalid, Valid}
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromwell.backend.BackendInitializationData
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.standard.{StandardExpressionFunctions, StandardInitializationActor, StandardInitializationActorParams}
import cromwell.backend.wfs.WorkflowPathBuilder
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._

import scala.concurrent.Future

class SharedFileSystemInitializationActor(standardParams: StandardInitializationActorParams)
  extends StandardInitializationActor(standardParams) {

  private implicit val system = context.system

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

  override lazy val pathBuilders: Future[List[PathBuilder]] =
    gcsPathBuilderFactory.toList.traverse(_.withOptions(workflowDescriptor.workflowOptions)).map(_ ++ Option(DefaultPathBuilder))

  override lazy val workflowPaths: Future[WorkflowPaths] = pathBuilders map {
    WorkflowPathBuilder.workflowPaths(configurationDescriptor, workflowDescriptor, _)
  }

  override lazy val expressionFunctions: Class[_ <: StandardExpressionFunctions] =
    classOf[SharedFileSystemExpressionFunctions]

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    initializationData map { data =>
      publishWorkflowRoot(data.workflowPaths.workflowRoot.pathAsString)
      data.workflowPaths.workflowRoot.createPermissionedDirectories()
      Option(data)
    }
  }
}
