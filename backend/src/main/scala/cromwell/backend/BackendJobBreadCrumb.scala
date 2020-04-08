package cromwell.backend

import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import cromwell.core.{JobKey, WorkflowId}
import wom.callable.Callable

case class BackendJobBreadCrumb(callable: Callable, id: WorkflowId, jobKey: JobKey) {
  def toPath(root: Path): Path = {
    val workflowPart = root.resolve(callable.name).resolve(id.toString)
    JobPaths.callPathBuilder(workflowPart, jobKey, isCallCacheCopyAttempt = false)
  }
}
