package cromwell.database.sql.joins

sealed trait MetadataJobQueryValue

/**
  * Applies a filter to retrieve only metadata keys matching this job
  */
case class CallQuery(callFqn: String, jobIndex: Option[Int], jobAttempt: Option[Int]) extends MetadataJobQueryValue

/**
  * Used to filter metadata keys at the workflow level only
  */
case object WorkflowQuery extends MetadataJobQueryValue

/**
  * Applies no filter on metadata job key. Both workflow level and job level keys are retrieved, for any job
  */
case object CallOrWorkflowQuery extends MetadataJobQueryValue
