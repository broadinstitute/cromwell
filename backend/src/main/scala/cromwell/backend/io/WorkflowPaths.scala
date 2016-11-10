package cromwell.backend.io

import java.nio.file.{FileSystem, FileSystems, Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.PathFactory
import net.ceedubs.ficus.Ficus._

class WorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, config: Config, val fileSystems: List[FileSystem] = List(FileSystems.getDefault)) extends PathFactory {
  val executionRoot = Paths.get(config.as[Option[String]]("root").getOrElse("cromwell-executions")).toAbsolutePath
  val DockerRoot = Paths.get(config.as[Option[String]]("dockerRoot").getOrElse("/root")).toAbsolutePath

  private def workflowPathBuilder(root: Path) = {
    root.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
        .resolve(workflowDescriptor.id.toString)
  }

  lazy val workflowRoot = workflowPathBuilder(executionRoot)
  lazy val dockerWorkflowRoot = workflowPathBuilder(DockerRoot)

  def toJobPaths(jobKey: BackendJobDescriptorKey) = new JobPaths(workflowDescriptor, config, jobKey)
}
