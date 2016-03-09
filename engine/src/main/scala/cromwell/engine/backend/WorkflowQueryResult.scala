package cromwell.engine.backend

import org.joda.time.DateTime

/** Corresponds to a single workflow returned within a `WorkflowQueryResponse`. */
case class WorkflowQueryResult(id: String, name: String, status: String, start: DateTime, end: Option[DateTime])
