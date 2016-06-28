package cromwell.backend.io

import java.nio.file.{Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.{BackendInitializationData, BackendWorkflowDescriptor}

object WorkflowPaths{
  val DockerRoot = Paths.get("/root")
}

class WorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, config: Config, initializationData: Option[BackendInitializationData]) {
  val executionRoot = Paths.get(config.getString("root"))

  private def workflowPathBuilder(root: Path) = {
    root.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
        .resolve(workflowDescriptor.id.toString)
  }

  lazy val workflowRoot = workflowPathBuilder(executionRoot)
  lazy val dockerWorkflowRoot = workflowPathBuilder(WorkflowPaths.DockerRoot)
}
