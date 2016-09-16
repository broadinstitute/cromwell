package cromwell.engine

import java.nio.file.{FileSystem, FileSystems}

import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.engine.backend.EnhancedWorkflowOptions._
import cromwell.filesystems.gcs.{GcsFileSystem, GcsFileSystemProvider, GoogleConfiguration}
import lenthall.config.ScalaConfig._
import lenthall.exception.MessageAggregation

import scala.concurrent.ExecutionContext

object EngineFilesystems {

  private val config = ConfigFactory.load
  private val googleConf: GoogleConfiguration = GoogleConfiguration(config)
  private val googleAuthMode = config.getStringOption("engine.filesystems.gcs.auth") map { confMode =>
    googleConf.auth(confMode) match {
      case scalaz.Success(mode) => mode
      case scalaz.Failure(errors) => throw new RuntimeException() with MessageAggregation {
        override def exceptionContext: String = s"Failed to create authentication mode for $confMode"
        override def errorMessages: Traversable[String] = errors.list.toList
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
