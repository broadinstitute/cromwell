package cromwell.backend.io

import java.nio.file.{Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.BackendWorkflowDescriptor

object WorkflowPaths{
  val DockerRoot = Paths.get("/root")
}

class WorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, config: Config) {
  val executionRoot = Paths.get(config.getString("root"))

  private def workflowPathBuilder(root: Path) = {
    root.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
        .resolve(workflowDescriptor.id.toString)
  }

  lazy val workflowRoot = workflowPathBuilder(executionRoot)
  lazy val dockerWorkflowRoot = workflowPathBuilder(WorkflowPaths.DockerRoot)
}
