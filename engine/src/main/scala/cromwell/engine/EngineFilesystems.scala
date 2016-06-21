package cromwell.engine

import java.nio.file.{FileSystem, FileSystems}

import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.engine.backend.EnhancedWorkflowOptions._
import cromwell.filesystems.gcs.{GcsFileSystemProvider, GoogleConfiguration, GcsFileSystem}
import lenthall.config.ScalaConfig._

import scala.util.Try

object EngineFilesystems {

  private val config = ConfigFactory.load

  def filesystemsForWorkflow(workflowOptions: WorkflowOptions): List[FileSystem] = {
    def gcsFileSystem: Option[GcsFileSystem] = {
      for {
        authModeString <- ConfigFactory.load.getStringOption("engine.filesystems.gcs.auth")
        authMode <- GoogleConfiguration(config).auth(authModeString).toOption
        fs <- Try(GcsFileSystem(GcsFileSystemProvider(
          authMode.buildStorage(workflowOptions.toGoogleAuthOptions, config)))).toOption
      } yield fs
    }

    List(gcsFileSystem, Option(FileSystems.getDefault)).flatten
  }
}
