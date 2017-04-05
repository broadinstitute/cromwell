package cromwell.api.model

import java.net.URL

/**
  * Represents information which we need to capture about a workflow sent to Cromwell.
  */
case class SubmittedWorkflow(id: WorkflowId, cromwellServer: URL, workflow: WorkflowSubmission)
