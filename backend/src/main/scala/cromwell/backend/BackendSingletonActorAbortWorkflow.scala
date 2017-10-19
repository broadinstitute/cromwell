package cromwell.backend

import cromwell.core.WorkflowId

case class BackendSingletonActorAbortWorkflow(workflowId: WorkflowId)
case object JobAborted
