package cromwell.backend.io

import java.nio.file.{FileSystems, FileSystem, Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.PathFactory

object WorkflowPaths{
  val DockerRoot = Paths.get("/root")
}

class WorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, config: Config, val fileSystems: List[FileSystem] = List(FileSystems.getDefault)) extends PathFactory {
  val executionRoot = Paths.get(config.getString("root"))

  private def workflowPathBuilder(root: Path) = {
    root.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
        .resolve(workflowDescriptor.id.toString)
  }

  lazy val workflowRoot = workflowPathBuilder(executionRoot)
  lazy val dockerWorkflowRoot = workflowPathBuilder(WorkflowPaths.DockerRoot)

  def toJobPaths(jobKey: BackendJobDescriptorKey) = new JobPaths(workflowDescriptor, config, jobKey)
}
