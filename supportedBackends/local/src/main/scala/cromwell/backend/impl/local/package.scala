package cromwell.backend.impl

import cromwell.backend.BackendInitializationData
import cromwell.backend.io.WorkflowPaths

package object local {
  case class LocalBackendInitializationData(workflowPaths: WorkflowPaths) extends BackendInitializationData
}
