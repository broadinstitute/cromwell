package cromwell.database.slick.tables

import java.sql.{Clob, Timestamp}

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

    def metadataValue = column[Option[Clob]]("METADATA_VALUE")

    def metadataValueType = column[Option[String]]("METADATA_VALUE_TYPE", O.Length(10))

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

  val metadataEntries = TableQuery[MetadataEntries]

  val metadataEntryIdsAutoInc = metadataEntries returning metadataEntries.map(_.metadataEntryId)

  val metadataEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val metadataEntryExistsForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).exists
  )

  val metadataEntriesForWorkflowExecutionUuidAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String], metadataKey: Rep[String]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.metadataKey === metadataKey
      if metadataEntry.callFullyQualifiedName.isEmpty
      if metadataEntry.jobIndex.isEmpty
      if metadataEntry.jobAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val metadataEntriesForJobKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFullyQualifiedName: Rep[String], jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Option[Int]]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
      if hasSameIndex(metadataEntry, jobIndex)
      if hasSameAttempt(metadataEntry, jobAttempt)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val metadataEntriesForJobKeyAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String], metadataKey: Rep[String], callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]], jobAttempt: Rep[Option[Int]]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.metadataKey === metadataKey
      if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
      if hasSameIndex(metadataEntry, jobIndex)
      if hasSameAttempt(metadataEntry, jobAttempt)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val metadataEntriesForIdGreaterThanOrEqual = Compiled(
    (metadataEntryId: Rep[Long], startMetadataKey: Rep[String], endMetadataKey: Rep[String], nameMetadataKey: Rep[String],
     statusMetadataKey: Rep[String], likeLabelMetadataKey: Rep[String]) => (for {
      metadataEntry <- metadataEntries
      if metadataEntry.metadataEntryId >= metadataEntryId
      if (metadataEntry.metadataKey === startMetadataKey || metadataEntry.metadataKey === endMetadataKey ||
        metadataEntry.metadataKey === nameMetadataKey || metadataEntry.metadataKey === statusMetadataKey ||
        metadataEntry.metadataKey.like(likeLabelMetadataKey)) &&
        (metadataEntry.callFullyQualifiedName.isEmpty && metadataEntry.jobIndex.isEmpty &&
          metadataEntry.jobAttempt.isEmpty)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  /**
    * Returns metadata entries that are "like" metadataKeys for she specified workflow.
    * If requireEmptyJobKey is true, only workflow level keys are returned, otherwise both workflow and call level
    * keys are returned.
    */
  def metadataEntriesLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String],
                                      requireEmptyJobKey: Boolean) = {
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeys)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  }

  /**
    * Returns metadata entries that are "like" metadataKeys for she specified call.
    * If jobAttempt has no value, all metadata keys for all attempts are returned.
    */
  def metadataEntriesLikeMetadataKeysWithJob(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String],
                                             callFqn: String, jobIndex: Option[Int], jobAttempt: Option[Int]) = {
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeys)
      if metadataEntry.callFullyQualifiedName === callFqn
      if hasSameIndex(metadataEntry, jobIndex)
      // Assume that every metadata entry for a call should have a non null attempt value
      // Because of that, if the jobAttempt paramater is Some(_), make sure it matches, otherwise take all entries
      // regardless of the attempt
      if (metadataEntry.jobAttempt === jobAttempt) || jobAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  }

  /**
    * Returns metadata entries that are NOT "like" metadataKeys for she specified workflow.
    * If requireEmptyJobKey is true, only workflow level keys are returned, otherwise both workflow and call level
    * keys are returned.
    */
  def metadataEntriesNotLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String],
                                         requireEmptyJobKey: Boolean) = {
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if !metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeys)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  }

  /**
    * Returns metadata entries that are NOT "like" metadataKeys for she specified call.
    * If jobIndex (resp. jobAttempt) has no value, all metadata keys for all indices (resp. attempt)
    * are returned.
    */
  def metadataEntriesNotLikeMetadataKeysWithJob(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String],
                                                callFqn: String, jobIndex: Option[Int], jobAttempt: Option[Int]) = {
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if !metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeys)
      if metadataEntry.callFullyQualifiedName === callFqn
      if hasSameIndex(metadataEntry, jobIndex)
        // Assume that every metadata entry for a call should have a non null attempt value
        // Because of that, if the jobAttempt paramater is Some(_), make sure it matches, otherwise take all entries
        // regardless of the attempt
      if (metadataEntry.jobAttempt === jobAttempt) || jobAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  }

  private[this] def metadataEntryHasMetadataKeysLike(metadataEntry: MetadataEntries,
                                                     metadataKeys: NonEmptyList[String]): Rep[Boolean] = {
    metadataKeys.toList.map(metadataEntry.metadataKey like _).reduce(_ || _)
  }

  private[this] def hasSameIndex(metadataEntry: MetadataEntries, jobIndex: Rep[Option[Int]]) = {
    (metadataEntry.jobIndex.isEmpty && jobIndex.isEmpty) || (metadataEntry.jobIndex === jobIndex)
  }

  private[this] def hasSameAttempt(metadataEntry: MetadataEntries, jobAttempt: Rep[Option[Int]]) = {
    (metadataEntry.jobAttempt.isEmpty && jobAttempt.isEmpty) || (metadataEntry.jobAttempt === jobAttempt)
  }

  private[this] def metadataEntryHasEmptyJobKey(metadataEntry: MetadataEntries,
                                                requireEmptyJobKey: Rep[Boolean]): Rep[Boolean] = {
    !requireEmptyJobKey ||
      (metadataEntry.callFullyQualifiedName.isEmpty &&
        metadataEntry.jobIndex.isEmpty &&
        metadataEntry.jobAttempt.isEmpty)
  }
}
