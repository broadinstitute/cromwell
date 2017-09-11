package cromwell

import cromwell.core.{JobKey, WorkflowId}
import cromwell.core.CromwellGraphNode._

package object jobstore {
  implicit class EnhancedJobKey(val jobKey: JobKey) extends AnyVal {
    def toJobStoreKey(workflowId: WorkflowId): JobStoreKey = JobStoreKey(workflowId, jobKey.scope.fullyQualifiedName, jobKey.index, jobKey.attempt)
  }
}
