package cromwell.backend.io

import java.nio.file.Paths

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.PathBuilder

object WorkflowPathsWithDocker {
  val DockerRoot = Paths.get("/root")
}

class WorkflowPathsWithDocker(val workflowDescriptor: BackendWorkflowDescriptor, val config: Config, val pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {
  val dockerWorkflowRoot = workflowPathBuilder(WorkflowPathsWithDocker.DockerRoot)
  override def toJobPaths(jobKey: BackendJobDescriptorKey): JobPaths = new JobPathsWithDocker(jobKey, workflowDescriptor, config, pathBuilders)
}