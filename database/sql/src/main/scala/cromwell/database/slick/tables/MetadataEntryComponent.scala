package cromwell.database.slick.tables

import java.sql.Timestamp

import cats.data.NonEmptyList
import cromwell.database.sql.tables.MetadataEntry

trait MetadataEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class MetadataEntries(tag: Tag) extends Table[MetadataEntry](tag, "METADATA_ENTRY") {
    /*
    DO NOT COPY/PASTE THIS CLASS!

    This class was unable to be standardized due to its size. The time to rewrite the rows/indexes was too large.

    http://stackoverflow.com/questions/2086105/why-does-it-take-so-long-to-rename-a-column-in-mysql

    Therefore, the scala database and sql database names DO NOT match, for now.

    Copy/paste another standardized table, and not this one.
     */

    def metadataEntryId = column[Long]("METADATA_JOURNAL_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID") // TODO: rename column via libquibase

    def callFullyQualifiedName = column[Option[String]]("CALL_FQN") // TODO: rename column via libquibase

    def jobIndex = column[Option[Int]]("JOB_SCATTER_INDEX") // TODO: rename column via libquibase

    def jobAttempt = column[Option[Int]]("JOB_RETRY_ATTEMPT") // TODO: rename column via libquibase

    def metadataKey = column[String]("METADATA_KEY")

    def metadataValue = column[Option[String]]("METADATA_VALUE")

    def metadataValueType = column[Option[String]]("METADATA_VALUE_TYPE")

    def metadataTimestamp = column[Timestamp]("METADATA_TIMESTAMP")

    override def * = (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt, metadataKey, metadataValue,
      metadataValueType, metadataTimestamp, metadataEntryId.?) <> (MetadataEntry.tupled, MetadataEntry.unapply)

    // TODO: rename index via libquibase
    def ixMetadataEntryWeu = index("METADATA_WORKFLOW_IDX", workflowExecutionUuid, unique = false)

    // TODO: rename index via libquibase
    def ixMetadataEntryWeuCfqnJiJa = index("METADATA_JOB_IDX",
      (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt), unique = false)

    // TODO: rename index via libquibase, and change column order from WEU_[MK]_CFQN_JI_JA to WEU_CFQN_JI_JA_[MK]
    def ixMetadataEntryWeuCfqnJiJaMk = index("METADATA_JOB_AND_KEY_IDX",
      (workflowExecutionUuid, metadataKey, callFullyQualifiedName, jobIndex, jobAttempt), unique = false)
  }

  protected val metadataEntries = TableQuery[MetadataEntries]

  val metadataEntryIdsAutoInc = metadataEntries returning metadataEntries.map(_.metadataEntryId)

  val metadataEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry
  )

  val metadataEntryExistsForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).exists
  )

  val metadataEntriesForWorkflowExecutionUuidAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String], metadataKey: Rep[String]) => for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.metadataKey === metadataKey
      if metadataEntry.callFullyQualifiedName.isEmpty
      if metadataEntry.jobIndex.isEmpty
      if metadataEntry.jobAttempt.isEmpty
    } yield metadataEntry
  )

  val metadataEntriesForJobKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFullyQualifiedName: Rep[String], jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Int]) => for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
      if (metadataEntry.jobIndex === jobIndex) ||
        (metadataEntry.jobIndex.isEmpty && jobIndex.isEmpty)
      if metadataEntry.jobAttempt === jobAttempt
    } yield metadataEntry
  )

  val metadataEntriesForJobKeyAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String], metadataKey: Rep[String], callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]], jobAttempt: Rep[Int]) => for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.metadataKey === metadataKey
      if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
      if (metadataEntry.jobIndex === jobIndex) ||
        (metadataEntry.jobIndex.isEmpty && jobIndex.isEmpty)
      if metadataEntry.jobAttempt === jobAttempt
    } yield metadataEntry
  )

  val metadataEntriesForIdGreaterThanOrEqual = Compiled(
    (metadataEntryId: Rep[Long], metadataKey1: Rep[String], metadataKey2: Rep[String], metadataKey3: Rep[String],
     metadataKey4: Rep[String]) => for {
      metadataEntry <- metadataEntries
      if metadataEntry.metadataEntryId >= metadataEntryId
      if (metadataEntry.metadataKey === metadataKey1 || metadataEntry.metadataKey === metadataKey2 ||
        metadataEntry.metadataKey === metadataKey3 || metadataEntry.metadataKey === metadataKey4) &&
        (metadataEntry.callFullyQualifiedName.isEmpty && metadataEntry.jobIndex.isEmpty &&
          metadataEntry.jobAttempt.isEmpty)
    } yield metadataEntry
  )

  def metadataEntriesLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String],
                                      requireEmptyJobKey: Boolean) = {
    for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeys)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry
  }

  def metadataEntriesNotLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String],
                                         requireEmptyJobKey: Boolean) = {
    for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if !metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeys)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry
  }

  private[this] def metadataEntryHasMetadataKeysLike(metadataEntry: MetadataEntries,
                                                     metadataKeys: NonEmptyList[String]): Rep[Boolean] = {
    metadataKeys.toList.map(metadataEntry.metadataKey like _).reduce(_ || _)
  }

  private[this] def metadataEntryHasEmptyJobKey(metadataEntry: MetadataEntries,
                                                requireEmptyJobKey: Rep[Boolean]): Rep[Boolean] = {
    !requireEmptyJobKey ||
      (metadataEntry.callFullyQualifiedName.isEmpty &&
        metadataEntry.jobIndex.isEmpty &&
        metadataEntry.jobAttempt.isEmpty)
  }
}
