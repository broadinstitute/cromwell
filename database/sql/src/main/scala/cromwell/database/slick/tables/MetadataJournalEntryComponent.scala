package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.MetadataJournalEntry
import javax.sql.rowset.serial.SerialClob

trait MetadataJournalEntryComponent {

  this: DriverComponent with WorkflowMetadataSummaryEntryComponent =>

  import driver.api._

  class MetadataJournalEntries(tag: Tag) extends Table[MetadataJournalEntry](tag, "METADATA_ENTRY") {

    def metadataJournalId = column[Long]("METADATA_JOURNAL_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255))

    def callFqn = column[Option[String]]("CALL_FQN", O.Length(255))

    def jobScatterIndex = column[Option[Int]]("JOB_SCATTER_INDEX")

    def jobRetryAttempt = column[Option[Int]]("JOB_RETRY_ATTEMPT")

    def metadataKey = column[String]("METADATA_KEY", O.Length(255))

    def metadataValue = column[Option[SerialClob]]("METADATA_VALUE")

    def metadataValueType = column[Option[String]]("METADATA_VALUE_TYPE", O.Length(10))

    def metadataTimestamp = column[Timestamp]("METADATA_TIMESTAMP")

    def needsSummarization = column[Boolean]("NEEDS_SUMMARIZATION")

    override def * = (workflowExecutionUuid, callFqn, jobScatterIndex, jobRetryAttempt, metadataKey, metadataValue,
      metadataValueType, metadataTimestamp, needsSummarization, metadataJournalId.?) <> (MetadataJournalEntry.tupled, MetadataJournalEntry.unapply)


    def ixMetadataJournalEntryWeu = index("IX_METADATA_JOURNAL_ENTRY_WEU", workflowExecutionUuid, unique = false)
    def ixMetadataJournalEntryNs = index("IX_METADATA_JOURNAL_ENTRY_NS", needsSummarization, unique = false)
  }

  val metadataJournalEntries = TableQuery[MetadataJournalEntries]

  val metadataJournalEntryIdsAutoInc = metadataJournalEntries returning metadataJournalEntries.map(_.metadataJournalId)

  val metadataJournalEntriesExists = Compiled(metadataJournalEntries.take(1).exists)

  val metadataJournalEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).sortBy(_.metadataTimestamp)
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
      if metadataEntry.callFqn.isEmpty
      if metadataEntry.jobScatterIndex.isEmpty
      if metadataEntry.jobRetryAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val metadataJournalEntriesForJobKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFullyQualifiedName: Rep[String], jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Option[Int]]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.callFqn === callFullyQualifiedName
      if hasSameIndex(metadataEntry, jobIndex)
      if hasSameAttempt(metadataEntry, jobAttempt)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val metadataJournalEntriesForJobKeyAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String], metadataKey: Rep[String], callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]], jobAttempt: Rep[Option[Int]]) => (for {
      metadataEntry <- metadataJournalEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntry.metadataKey === metadataKey
      if metadataEntry.callFqn === callFullyQualifiedName
      if hasSameIndex(metadataEntry, jobIndex)
      if hasSameAttempt(metadataEntry, jobAttempt)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  def entriesInNeedOfSummarization(limit: Int) = Compiled(
    () => (for {
      entry <- metadataJournalEntries
      if entry.needsSummarization
    } yield entry).sortBy(_.metadataJournalId).take(limit)
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
    } yield metadataEntry).sortBy(_.metadataTimestamp)
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
      if metadataEntry.callFqn === callFqn
      if hasSameIndex(metadataEntry, jobIndex)
      // Assume that every metadata entry for a call should have a non null attempt value
      // Because of that, if the jobAttempt parameter is Some(_), make sure it matches, otherwise take all entries
      // regardless of the attempt
      if (metadataEntry.jobRetryAttempt === jobAttempt) || jobAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataTimestamp)
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
    (metadataEntry.jobScatterIndex.isEmpty && jobIndex.isEmpty) || (metadataEntry.jobScatterIndex === jobIndex)
  }

  private[this] def hasSameAttempt(metadataEntry: MetadataJournalEntries, jobAttempt: Rep[Option[Int]]) = {
    (metadataEntry.jobRetryAttempt.isEmpty && jobAttempt.isEmpty) || (metadataEntry.jobRetryAttempt === jobAttempt)
  }

  private[this] def metadataEntryHasEmptyJobKey(metadataEntry: MetadataJournalEntries,
                                                requireEmptyJobKey: Rep[Boolean]): Rep[Boolean] = {
    !requireEmptyJobKey ||
      (metadataEntry.callFqn.isEmpty &&
        metadataEntry.jobScatterIndex.isEmpty &&
        metadataEntry.jobRetryAttempt.isEmpty)
  }
}
