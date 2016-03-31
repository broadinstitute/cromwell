package cromwell.engine.analysis

import cromwell.engine.backend.{WorkflowDescriptor, Backend}
import wdl4s.Scope

// This can be enhanced later...
object BackendSelector {
  def selectBackend(workflowDescriptor: WorkflowDescriptor, scope: Scope): Backend =
    // TODO: Probably want to enhance this somewhat:
    workflowDescriptor.defaultBackend
}
