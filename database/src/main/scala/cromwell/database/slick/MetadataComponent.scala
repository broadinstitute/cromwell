package cromwell.database.slick


import java.sql.Timestamp

import cromwell.database.obj.{Metadatum, WorkflowMetadataKeys}
import scalaz._

trait MetadataComponent {

  this: DriverComponent =>

  import driver.api._

  class Metadata(tag: Tag) extends Table[Metadatum](tag, "METADATA_JOURNAL") {
    def metadataId = column[Long]("METADATA_JOURNAL_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")
    def key = column[String]("METADATA_KEY")
    def callFqn = column[Option[String]]("METADATA_CALL_FQN")
    def index = column[Option[Int]]("METADATA_CALL_INDEX")
    def attempt = column[Option[Int]]("METADATA_CALL_ATTEMPT")
    def value = column[Option[String]]("METADATA_VALUE")
    def valueType = column[Option[String]]("METADATA_VALUE_TYPE")
    def timestamp = column[Timestamp]("METADATA_TIMESTAMP")

    override def * = (workflowExecutionUuid, key, callFqn, index, attempt, value, valueType, timestamp, metadataId.?) <>
      (Metadatum.tupled, Metadatum.unapply)

    def workflowIndex = index("METADATA_WORKFLOW_IDX", workflowExecutionUuid, unique = false)

    def jobIndex = index("METADATA_JOB_IDX", (workflowExecutionUuid, callFqn, index, attempt), unique = false)

    def jobAndKeyIndex = index("METADATA_JOB_AND_KEY_IDX", (workflowExecutionUuid, key, callFqn, index, attempt), unique = false)
  }

  protected val metadata = TableQuery[Metadata]

  val metadataByWorkflowUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      metadatum <- metadata
      if metadatum.workflowExecutionUuid === workflowExecutionUuid
    } yield metadatum)

  val metadataByWorkflowUuidAndKey = Compiled(
    (workflowExecutionUuid: Rep[String], key: Rep[String]) => for {
      metadatum <- metadata
      if metadatum.workflowExecutionUuid === workflowExecutionUuid
      if metadatum.key === key
      if metadatum.callFqn.isEmpty
      if metadatum.index.isEmpty
      if metadatum.attempt.isEmpty
    } yield metadatum)

  val metadataByWorkflowUuidAndCallFqnAndIndexAndAttempt = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], index: Rep[Option[Int]], attempt: Rep[Int]) => for {
      metadatum <- metadata
      if metadatum.workflowExecutionUuid === workflowExecutionUuid
      if metadatum.callFqn === callFqn
      if (metadatum.index === index) || (metadatum.index.isEmpty && index.isEmpty)
      if metadatum.attempt === attempt
    } yield metadatum)

  val metadataByWorkflowUuidAndKeyAndCallFqnAndIndexAndAttempt = Compiled(
    (workflowExecutionUuid: Rep[String], key: Rep[String], callFqn: Rep[String], index: Rep[Option[Int]], attempt: Rep[Int]) => for {
      metadatum <- metadata
      if metadatum.workflowExecutionUuid === workflowExecutionUuid
      if metadatum.key === key
      if metadatum.callFqn === callFqn
      if (metadatum.index === index) || (metadatum.index.isEmpty && index.isEmpty)
      if metadatum.attempt === attempt
    } yield metadatum)

  val metadataWithIdAndTimestampGreaterThanOrEqual = Compiled(
    (id: Rep[Long], timestamp: Rep[Option[Timestamp]]) => for {
      m <- metadata
      if m.metadataId >= id && (m.timestamp >= timestamp || timestamp.isEmpty)
      if (m.key === WorkflowMetadataKeys.StartTime ||
        m.key === WorkflowMetadataKeys.EndTime ||
        m.key === WorkflowMetadataKeys.Name ||
        m.key === WorkflowMetadataKeys.Status) &&
        (m.callFqn.isEmpty && m.index.isEmpty && m.attempt.isEmpty)
    } yield m
  )

  val metadataAutoInc = metadata returning metadata.map(_.metadataId)

  def queryMetadataMatchingAnyWildcardKeys(workflowUuid: String, keys: NonEmptyList[String], requireEmptyJobKey: Boolean) = {
    def hasEmptyJobKey(m: Metadata): Rep[Boolean] = m.callFqn.isEmpty && m.index.isEmpty && m.attempt.isEmpty
    val repRequireEmptyJobKey: Rep[Boolean] = requireEmptyJobKey

    for {
      m <- metadata
      if m.workflowExecutionUuid === workflowUuid
      if keys.list map { m.key like _ } reduce (_ || _)
      if !repRequireEmptyJobKey || hasEmptyJobKey(m)
    } yield m
  }
}
