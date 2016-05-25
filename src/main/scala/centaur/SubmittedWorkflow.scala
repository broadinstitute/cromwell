package centaur

import java.net.URL
import java.util.UUID

import centaur.test.workflow.Workflow


/**
  * Represents information which we need to capture about a workflow sent to Cromwell.
  */
case class SubmittedWorkflow(id: UUID, cromwellServer: URL, workflow: Workflow)

