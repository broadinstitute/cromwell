package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder}

object WorkflowPathsWithDocker {
  val DockerRoot: Path = DefaultPathBuilder.get("/cromwell-executions")
}

class WorkflowPathsWithDocker(val workflowDescriptor: BackendWorkflowDescriptor, val config: Config, val pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {
  val dockerWorkflowRoot: Path = workflowPathBuilder(WorkflowPathsWithDocker.DockerRoot)

  override def toJobPaths(jobKey: BackendJobDescriptorKey, jobWorkflowDescriptor: BackendWorkflowDescriptor): JobPathsWithDocker = new JobPathsWithDocker(this, jobKey)
}
