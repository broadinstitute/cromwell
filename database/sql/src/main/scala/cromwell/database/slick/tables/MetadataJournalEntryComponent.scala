package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.{MetadataEntry, MetadataJournalEntry}
import javax.sql.rowset.serial.SerialClob

trait MetadataJournalEntryComponent {

  this: DriverComponent with WorkflowMetadataSummaryEntryComponent =>

  import driver.api._

  class MetadataJournalEntries(tag: Tag) extends Table[MetadataJournalEntry](tag, "METADATA_ENTRY") {
    /*
    DO NOT COPY/PASTE THIS CLASS!

    This class was unable to be standardized due to its size. The time to rewrite the rows/indexes was too large.

    http://stackoverflow.com/questions/2086105/why-does-it-take-so-long-to-rename-a-column-in-mysql

    Therefore, the scala database and sql database names DO NOT match, for now.

    Copy/paste another standardized table, and not this one.
     */

    def metadataJournalEntryId = column[Long]("METADATA_JOURNAL_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255)) // TODO: rename column via liquibase

    def callFullyQualifiedName = column[Option[String]]("CALL_FQN", O.Length(255)) // TODO: rename column via liquibase

    def jobIndex = column[Option[Int]]("JOB_SCATTER_INDEX") // TODO: rename column via liquibase

    def jobAttempt = column[Option[Int]]("JOB_RETRY_ATTEMPT") // TODO: rename column via liquibase

    def metadataKey = column[String]("METADATA_KEY", O.Length(255))

    def metadataValue = column[Option[SerialClob]]("METADATA_VALUE")

    def metadataValueType = column[Option[String]]("METADATA_VALUE_TYPE", O.Length(10))

    def metadataGenerationTimestamp = column[Timestamp]("METADATA_GENERATION_TIMESTAMP")

    def metadataWriteTimestamp = column[Timestamp]("METADATA_WRITE_TIMESTAMP")

    def summarizationTimestamp = column[Option[Timestamp]]("SUMMARIZATION_TIMESTAMP")

    override def * = (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt, metadataKey, metadataValue,
      metadataValueType, metadataGenerationTimestamp, metadataWriteTimestamp, summarizationTimestamp, metadataJournalEntryId.?) <> (MetadataJournalEntry.tupled, MetadataJournalEntry.unapply)

    // TODO: rename index via liquibase
    def ixMetadataEntryWeu = index("METADATA_WORKFLOW_IDX", workflowExecutionUuid, unique = false)
  }

  val metadataJournalEntries = TableQuery[MetadataJournalEntries]

  val metadataJournalEntryIdsAutoInc = metadataJournalEntries returning metadataJournalEntries.map(_.metadataJournalEntryId)

  val metadataJournalEntriesExists = Compiled(metadataJournalEntries.take(1).exists)

  val metadataJournalEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).sortBy(_.metadataGenerationTimestamp)
  )

  val metadataJournalEntriesWithoutLabelsForRootWorkflowId = Compiled(
    (rootWorkflowId: Rep[String]) => {
      val targetWorkflowIds = for {
        summary <- workflowMetadataSummaryEntries
        // Uses `IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU`, `UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU`
        if summary.rootWorkflowExecutionUuid === rootWorkflowId || summary.workflowExecutionUuid === rootWorkflowId
      } yield summary.workflowExecutionUuid

      for {
        metadata <- metadataJournalEntries
        if metadata.workflowExecutionUuid in targetWorkflowIds // Uses `METADATA_WORKFLOW_IDX`
        if !(metadata.metadataKey like "labels:%")
      } yield metadata
    }
  )

  val metadataJournalEntryExistsForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).exists
  )

  def metadataJournalEntryExistsForWorkflowExecutionUuid(workflowId: Rep[String], key: Rep[String]): Rep[Boolean] = {
    metadataJournalEntries.filter( metadataEntry =>
      metadataEntry.workflowExecutionUuid === workflowId &&
      metadataEntry.metadataKey === key &&
      metadataEntry.metadataValue.isDefined
    ).exists
  }

  val metadataJournalEntriesForWorkflowExecutionUuidAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String], metadataKey: Rep[String]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.metadataKey === metadataKey
      if metadataEntry.callFullyQualifiedName.isEmpty
      if metadataEntry.jobIndex.isEmpty
      if metadataEntry.jobAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataGenerationTimestamp)
  )

  val metadataJournalEntriesForJobKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFullyQualifiedName: Rep[String], jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Option[Int]]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
      if hasSameIndex(metadataEntry, jobIndex)
      if hasSameAttempt(metadataEntry, jobAttempt)
    } yield metadataEntry).sortBy(_.metadataGenerationTimestamp)
  )

  val metadataJournalEntriesForJobKeyAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String], metadataKey: Rep[String], callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]], jobAttempt: Rep[Option[Int]]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.metadataKey === metadataKey
      if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
      if hasSameIndex(metadataEntry, jobIndex)
      if hasSameAttempt(metadataEntry, jobAttempt)
    } yield metadataEntry).sortBy(_.metadataGenerationTimestamp)
  )

  val metadataJournalEntriesForIdRange = Compiled(
    (minMetadataEntryId: Rep[Long], maxMetadataEntryId: Rep[Long]) => {
      for {
        metadataEntry <- metadataJournalEntries
        if metadataEntry.metadataJournalEntryId >= minMetadataEntryId
        if metadataEntry.metadataJournalEntryId <= maxMetadataEntryId
      } yield metadataEntry
    }
  )

  /**
    * Returns metadata entries that are "like" metadataKeys for the specified workflow.
    * If requireEmptyJobKey is true, only workflow level keys are returned, otherwise both workflow and call level
    * keys are returned.
    */
  def metadataJournalEntriesWithKeyConstraints(workflowExecutionUuid: String,
                                        metadataKeysToFilterFor: List[String],
                                        metadataKeysToFilterOut: List[String],
                                        requireEmptyJobKey: Boolean) = {
    (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataJournalEntryHasMetadataKeysLike(metadataEntry, metadataKeysToFilterFor, metadataKeysToFilterOut)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry).sortBy(_.metadataGenerationTimestamp)
  }

  /**
    * Returns metadata entries that are "like" metadataKeys for the specified call.
    * If jobAttempt has no value, all metadata keys for all attempts are returned.
    */
  def metadataJournalEntriesForJobWithKeyConstraints(workflowExecutionUuid: String,
                                              metadataKeysToFilterFor: List[String],
                                              metadataKeysToFilterOut: List[String],
                                              callFqn: String,
                                              jobIndex: Option[Int],
                                              jobAttempt: Option[Int]) = {
    (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataJournalEntryHasMetadataKeysLike(metadataEntry, metadataKeysToFilterFor, metadataKeysToFilterOut)
      if metadataEntry.callFullyQualifiedName === callFqn
      if hasSameIndex(metadataEntry, jobIndex)
      // Assume that every metadata entry for a call should have a non null attempt value
      // Because of that, if the jobAttempt parameter is Some(_), make sure it matches, otherwise take all entries
      // regardless of the attempt
      if (metadataEntry.jobAttempt === jobAttempt) || jobAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataGenerationTimestamp)
  }

  private[this] def metadataJournalEntryHasMetadataKeysLike(metadataEntry: MetadataJournalEntries,
                                                     metadataKeysToFilterFor: List[String],
                                                     metadataKeysToFilterOut: List[String]): Rep[Boolean] = {

    def containsKey(key: String): Rep[Boolean] = metadataEntry.metadataKey like key

    val positiveFilter: Option[Rep[Boolean]] = metadataKeysToFilterFor.map(containsKey).reduceOption(_ || _)
    val negativeFilter: Option[Rep[Boolean]] = metadataKeysToFilterOut.map(containsKey).reduceOption(_ || _)

    (positiveFilter, negativeFilter) match {
      case (Some(pf), Some(nf)) => pf && !nf
      case (Some(pf), None) => pf
      case (None, Some(nf)) => !nf

      // We should never get here, but there's no reason not to handle it:
      // ps: is there a better literal "true" in slick?
      case (None, None) => true: Rep[Boolean]
    }
  }

  private[this] def hasSameIndex(metadataEntry: MetadataJournalEntries, jobIndex: Rep[Option[Int]]) = {
    (metadataEntry.jobIndex.isEmpty && jobIndex.isEmpty) || (metadataEntry.jobIndex === jobIndex)
  }

  private[this] def hasSameAttempt(metadataEntry: MetadataJournalEntries, jobAttempt: Rep[Option[Int]]) = {
    (metadataEntry.jobAttempt.isEmpty && jobAttempt.isEmpty) || (metadataEntry.jobAttempt === jobAttempt)
  }

  private[this] def metadataEntryHasEmptyJobKey(metadataEntry: MetadataJournalEntries,
                                                requireEmptyJobKey: Rep[Boolean]): Rep[Boolean] = {
    !requireEmptyJobKey ||
      (metadataEntry.callFullyQualifiedName.isEmpty &&
        metadataEntry.jobIndex.isEmpty &&
        metadataEntry.jobAttempt.isEmpty)
  }
}
