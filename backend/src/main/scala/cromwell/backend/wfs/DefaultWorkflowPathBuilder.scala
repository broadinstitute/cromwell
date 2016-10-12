package cromwell.backend.wfs

import cromwell.core.path.DefaultPathBuilder


object DefaultWorkflowPathBuilder extends WorkflowPathBuilder {
  override def pathBuilderOption(params: WorkflowFileSystemProviderParams) = Option(DefaultPathBuilder)
}
