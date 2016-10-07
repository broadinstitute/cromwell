package cromwell.engine

import java.nio.file.{FileSystem, FileSystems}

import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.engine.backend.EnhancedWorkflowOptions._
import cromwell.filesystems.gcs.{GcsFileSystem, GcsFileSystemProvider, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext

object EngineFilesystems {

  private val config = ConfigFactory.load
  private val googleConf: GoogleConfiguration = GoogleConfiguration(config)
  private val googleAuthMode = config.as[Option[String]]("engine.filesystems.gcs.auth") map { confMode =>
    googleConf.auth(confMode) match {
      case Valid(mode) => mode
      case Invalid(errors) => throw new RuntimeException() with MessageAggregation {
        override def exceptionContext: String = s"Failed to create authentication mode for $confMode"
        override def errorMessages: Traversable[String] = errors.toList
      }
    }
  }

  def filesystemsForWorkflow(workflowOptions: WorkflowOptions)(implicit ec: ExecutionContext): List[FileSystem] = {
    def gcsFileSystem: Option[GcsFileSystem] = {
      googleAuthMode map { mode =>
        val storage = mode.buildStorage(workflowOptions.toGoogleAuthOptions, googleConf.applicationName)
        GcsFileSystem(GcsFileSystemProvider(storage))
      }
    }

    List(gcsFileSystem, Option(FileSystems.getDefault)).flatten
  }
}
