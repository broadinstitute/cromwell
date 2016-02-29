package centaur

import java.net.URL
import java.util.UUID

/**
  * Represents information which we need to capture about a workflow sent to Cromwell.
  */
case class Workflow(id: UUID, cromwellServer: URL)

