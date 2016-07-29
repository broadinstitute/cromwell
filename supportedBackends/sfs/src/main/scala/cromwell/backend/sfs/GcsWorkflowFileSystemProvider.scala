package cromwell.backend.sfs

import cromwell.backend.wfs.{WorkflowFileSystemProvider, WorkflowFileSystemProviderParams}
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions
import cromwell.filesystems.gcs.{GcsFileSystem, GcsFileSystemProvider, GoogleConfiguration}
import lenthall.config.ScalaConfig._
import wdl4s.ValidationException

import scala.util.Try

object GcsWorkflowFileSystemProvider extends WorkflowFileSystemProvider {
  override def fileSystemOption(params: WorkflowFileSystemProviderParams): Option[GcsFileSystem] = {
    params.fileSystemConfig.getStringOption("gcs.auth") map gcsFileSystem(params)
  }

  private def gcsFileSystem(params: WorkflowFileSystemProviderParams)(gcsAuthName: String): GcsFileSystem = {
    val workflowOptions = params.workflowOptions
    val globalConfig = params.globalConfig
    val googleConfig = GoogleConfiguration(globalConfig)
    val googleAuthModeValidation = googleConfig.auth(gcsAuthName)

    val gcsAuthMode = googleAuthModeValidation match {
      case scalaz.Success(googleAuthMode) => googleAuthMode
      case scalaz.Failure(errors) =>
        throw new ValidationException("Could not create gcs filesystem from configuration", errors)
    }

    val authOptions = new GoogleAuthOptions {
      override def get(key: String): Try[String] = workflowOptions.get(key)
    }

    val storage = gcsAuthMode.buildStorage(authOptions, googleConfig)
    GcsFileSystem(GcsFileSystemProvider(storage)(params.fileSystemExecutionContext))
  }
}
