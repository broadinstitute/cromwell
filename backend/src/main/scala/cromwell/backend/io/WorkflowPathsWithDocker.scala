package cromwell.backend.io

import java.nio.file.{Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.PathBuilder

object WorkflowPathsWithDocker {
  val DockerRoot: Path = Paths.get("/root")
}

class WorkflowPathsWithDocker(val workflowDescriptor: BackendWorkflowDescriptor, val config: Config, val pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {
  val dockerWorkflowRoot: Path = workflowPathBuilder(WorkflowPathsWithDocker.DockerRoot)

  override def toJobPaths(jobKey: BackendJobDescriptorKey,
                          jobWorkflowDescriptor: BackendWorkflowDescriptor): JobPathsWithDocker = {
    new JobPathsWithDocker(jobKey, jobWorkflowDescriptor, config, pathBuilders)
  }
}
