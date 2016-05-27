package cromwell.engine.backend

import java.time.OffsetDateTime

/** Corresponds to a single workflow returned within a `WorkflowQueryResponse`. */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class WorkflowQueryResult(id: String, name: String, status: String, start: OffsetDateTime, end: Option[OffsetDateTime])
