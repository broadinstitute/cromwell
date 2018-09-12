package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder}
import net.ceedubs.ficus.Ficus._

object WorkflowPathsWithDocker {
  val DockerRoot: Path = DefaultPathBuilder.get("/cromwell-executions")
}

final case class WorkflowPathsWithDocker(workflowDescriptor: BackendWorkflowDescriptor, config: Config, pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {
  val dockerRoot: Path = config.as[Option[String]]("dockerRoot") match {
    case Some(pathString) => DefaultPathBuilder.get(pathString)
    case _ => WorkflowPathsWithDocker.DockerRoot
  }
  val dockerWorkflowRoot: Path = workflowPathBuilder(dockerRoot)

  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): JobPathsWithDocker = {
    new JobPathsWithDocker(workflowPaths.asInstanceOf[WorkflowPathsWithDocker], jobKey)
  }

  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)
}
