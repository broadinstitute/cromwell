package cromwell.backend.wfs

import com.typesafe.config.Config
import cromwell.backend.io.{WorkflowPathsWithDocker, WorkflowPaths}
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilder

import scala.concurrent.ExecutionContext

object WorkflowPathBuilder {
  def workflowPaths(configurationDescriptor: BackendConfigurationDescriptor,
                    workflowDescriptor: BackendWorkflowDescriptor,
                    pathBuilders: List[PathBuilder]): WorkflowPaths = {
    new WorkflowPathsWithDocker(workflowDescriptor, configurationDescriptor.backendConfig, pathBuilders)
  }
}

final case class WorkflowFileSystemProviderParams(fileSystemConfig: Config, globalConfig: Config,
                                                  workflowOptions: WorkflowOptions,
                                                  fileSystemExecutionContext: ExecutionContext)

trait WorkflowPathBuilder {
  def pathBuilderOption(params: WorkflowFileSystemProviderParams): Option[PathBuilder]
}
