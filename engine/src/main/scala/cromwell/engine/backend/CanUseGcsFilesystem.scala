package cromwell.engine.backend

import java.nio.file.FileSystem

import cromwell.core.WorkflowOptions
import cromwell.engine.backend.EnhancedWorkflowOptions._
import cromwell.filesystems.gcs.{GcsFileSystemProvider, GoogleAuthMode, GoogleConfiguration}
import lenthall.config.ScalaConfig._
import lenthall.exception.MessageAggregation

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait CanUseGcsFilesystem { self: OldStyleBackend =>

  protected def gcsFilesystem(options: WorkflowOptions): Option[FileSystem] = {

    def gcsFilesystemMustHaveAuth: Option[String] = {
      backendConfig.getConfigOption("filesystems.gcs") map {
        _.getStringOr("auth", throw new RuntimeException(s"No auth in filesystems.gcs for backend '$name'.")) }
    }

    def gcsFilesystemAuthMustBeLegit(authName: String): GoogleAuthMode = {
      class ErrorMessageAggregatedException(override val exceptionContext: String,
                                            override val errorMessages: Traversable[String]) extends MessageAggregation

      GoogleConfiguration(globalConfig).auth(authName) valueOr { errs => throw new ErrorMessageAggregatedException(s"For backend '$name': ", errs.list) }
    }

    val gcsAuth: Option[GoogleAuthMode] = for {
      // Throw if there's a filesystems.gcs but no auth within that.
      authName <- gcsFilesystemMustHaveAuth
      // Throw if there's an auth named in this backend's `filesystems.gcs` which isn't in the global `google` list of auths.
      auth = gcsFilesystemAuthMustBeLegit(authName)
    } yield auth

    val maybeStorage = gcsAuth map { _.buildStorage(options.toGoogleAuthOptions, globalConfig) }
    maybeStorage map GcsFileSystemProvider.apply map { _.getFileSystem }
  }
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait MustHaveGcsFilesystem extends CanUseGcsFilesystem { self: OldStyleBackend =>
  override protected def gcsFilesystem(options: WorkflowOptions): Option[FileSystem] = {
    super.gcsFilesystem(options).orElse(
      throw new RuntimeException(s"Could not build GCS filesystem for backend '$name'"))
  }
}
