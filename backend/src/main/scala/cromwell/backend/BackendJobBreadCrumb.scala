package cromwell.backend

import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import cromwell.core.{JobKey, WorkflowId}
import wdl4s.wom.callable.WorkflowDefinition

case class BackendJobBreadCrumb(workflow: WorkflowDefinition, id: WorkflowId, jobKey: JobKey) {
  def toPath(root: Path): Path = {
    val workflowPart = root.resolve(workflow.name).resolve(id.toString)
    JobPaths.callPathBuilder(workflowPart, jobKey)
  }
}
