package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.Metadatum

import scalaz._

trait MetadataComponent {

  this: DriverComponent =>

  import driver.api._

  class Metadata(tag: Tag) extends Table[Metadatum](tag, "METADATA_JOURNAL") {
    def metadataId = column[Long]("METADATA_JOURNAL_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")
    def key = column[String]("METADATA_KEY")
    def callFqn = column[Option[String]]("CALL_FQN")
    def index = column[Option[Int]]("JOB_SCATTER_INDEX")
    def attempt = column[Option[Int]]("JOB_RETRY_ATTEMPT")
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

  val metadataWithIdGreaterThanOrEqual = Compiled(
    (id: Rep[Long], key1: Rep[String], key2: Rep[String], key3: Rep[String], key4: Rep[String]) => for {
      m <- metadata
      if m.metadataId >= id
      if (m.key === key1 || m.key === key2 || m.key === key3 || m.key === key4 ) &&
        (m.callFqn.isEmpty && m.index.isEmpty && m.attempt.isEmpty)
    } yield m
  )

  val workflowExists = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      m <- metadata
      if m.workflowExecutionUuid === workflowExecutionUuid
    } yield m).exists
  )

  val metadataAutoInc = metadata returning metadata.map(_.metadataId)

  def queryMetadataMatchingAnyWildcardKeys(workflowUuid: String, keys: NonEmptyList[String], requireEmptyJobKey: Boolean) = {
    val repRequireEmptyJobKey: Rep[Boolean] = requireEmptyJobKey

    for {
      m <- metadata
      if m.workflowExecutionUuid === workflowUuid
      if keys.list.toList map { m.key like _ } reduce (_ || _)
      if !repRequireEmptyJobKey || hasEmptyJobKey(m)
    } yield m
  }

  def queryMetadataNotMatchingAnyWildcardKeys(workflowUuid: String, keys: NonEmptyList[String],
                                              requireEmptyJobKey: Boolean) = {
    val repRequireEmptyJobKey: Rep[Boolean] = requireEmptyJobKey

    for {
      metadatum <- metadata
      if metadatum.workflowExecutionUuid === workflowUuid
      if !(keys.list.toList map { metadatum.key like _ } reduce (_ || _))
      if !repRequireEmptyJobKey || hasEmptyJobKey(metadatum)
    } yield metadatum
  }

  private[this] def hasEmptyJobKey(metadatum: Metadata): Rep[Boolean] = {
    metadatum.callFqn.isEmpty && metadatum.index.isEmpty && metadatum.attempt.isEmpty
  }

}
