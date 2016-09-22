package cromwell.backend.wfs

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilder
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext

object WorkflowPathBuilder {
  def workflowPaths(configurationDescriptor: BackendConfigurationDescriptor,
                    workflowDescriptor: BackendWorkflowDescriptor,
                    pathBuilders: List[PathBuilder],
                    fileSystemExecutionContext: ExecutionContext): WorkflowPaths = {
    val backendConfig = configurationDescriptor.backendConfig
    val fileSystemConfig = backendConfig.as[Option[Config]]("filesystems").getOrElse(ConfigFactory.empty())
    val globalConfig = configurationDescriptor.globalConfig
    val params = WorkflowFileSystemProviderParams(fileSystemConfig, globalConfig, workflowDescriptor.workflowOptions,
      fileSystemExecutionContext)
    new WorkflowPaths(workflowDescriptor, configurationDescriptor.backendConfig, pathBuilders)
  }
}

final case class WorkflowFileSystemProviderParams(fileSystemConfig: Config, globalConfig: Config,
                                                  workflowOptions: WorkflowOptions,
                                                  fileSystemExecutionContext: ExecutionContext)

trait WorkflowPathBuilder {
  def pathBuilderOption(params: WorkflowFileSystemProviderParams): Option[PathBuilder]
}
