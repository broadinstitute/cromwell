package cromwell.backend

import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import cromwell.core.{JobKey, WorkflowId}
import wdl4s.wdl.WdlWorkflow

case class BackendJobBreadCrumb(workflow: WdlWorkflow, id: WorkflowId, jobKey: JobKey) {
  def toPath(root: Path): Path = {
    val workflowPart = root.resolve(workflow.unqualifiedName).resolve(id.toString)
    JobPaths.callPathBuilder(workflowPart, jobKey)
  }
}
