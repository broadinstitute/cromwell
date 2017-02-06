package cromwell.backend

import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import cromwell.core.{JobKey, WorkflowId}
import wdl4s.Workflow

case class BackendJobBreadCrumb(workflow: Workflow, id: WorkflowId, jobKey: JobKey) {
  def toPath(root: Path): Path = {
    val workflowPart = root.resolve(workflow.unqualifiedName).resolve(id.toString)
    JobPaths.callPathBuilder(workflowPart, jobKey)
  }
}
