package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder}
import net.ceedubs.ficus.Ficus._

object WorkflowPathsWithDocker {
  val DefaultDockerRoot = "/cromwell-executions"
}

final case class WorkflowPathsWithDocker(workflowDescriptor: BackendWorkflowDescriptor, config: Config, pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {
  val dockerRoot: Path =
    DefaultPathBuilder.get(
      config.getOrElse[String]("dockerRoot", WorkflowPathsWithDocker.DefaultDockerRoot)
    )
  val dockerWorkflowRoot: Path = workflowPathBuilder(dockerRoot)

  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): JobPathsWithDocker = {
    new JobPathsWithDocker(workflowPaths.asInstanceOf[WorkflowPathsWithDocker], jobKey, isCallCacheCopyAttempt = false)
  }

  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)
}
