package cromwell.backend.sfs

import better.files._
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.BackendInitializationData
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.standard.{DefaultInitializationActorParams, StandardInitializationActor, StandardInitializationActorParams, StandardInitializationData}
import cromwell.backend.wfs.WorkflowPathBuilder
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilder, PathBuilderFactory}
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._

import scala.concurrent.Future
import scala.util.Try

/**
  * Initializes a shared file system actor factory and creates initialization data to pass to the execution actors.
  *
  * @param params Initialization parameters.
  */
class SharedFileSystemInitializationActor(params: DefaultInitializationActorParams)
  extends StandardInitializationActor {

  override val standardParams: StandardInitializationActorParams = params

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

  lazy val pathBuilders: List[PathBuilder] =
    pathBuilderFactories map { _.withOptions(workflowDescriptor.workflowOptions)(context.system) }

  val workflowPaths: WorkflowPaths =
    WorkflowPathBuilder.workflowPaths(configurationDescriptor, workflowDescriptor, pathBuilders)

  def initializationData: StandardInitializationData = {
    new StandardInitializationData(workflowPaths, runtimeAttributesBuilder)
  }

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.toString)
      File(workflowPaths.workflowRoot).createDirectories()
      Option(initializationData)
    })
  }
}
