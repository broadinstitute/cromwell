package cromwell.backend.sfs

import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.wfs.{WorkflowFileSystemProviderParams, WorkflowPathBuilder}
import cromwell.core.path.PathBuilder
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.config.ScalaConfig._
import wdl4s.ValidationException

object GcsWorkflowPathBuilder extends WorkflowPathBuilder {
  override def pathBuilderOption(params: WorkflowFileSystemProviderParams): Option[PathBuilder] = {
    params.fileSystemConfig.getStringOption("gcs.auth") map pathBuilder(params)
  }

  private def pathBuilder(params: WorkflowFileSystemProviderParams)(gcsAuthName: String): PathBuilder = {
    GoogleConfiguration(params.globalConfig).auth(gcsAuthName) match {
      case Valid(authMode) => GcsPathBuilderFactory(authMode).withOptions(params.workflowOptions)
      case Invalid(errors) => throw new ValidationException("Could not create gcs filesystem from configuration", errors)
    }
  }
}
