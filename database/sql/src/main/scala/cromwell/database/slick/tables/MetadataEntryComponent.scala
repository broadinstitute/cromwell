package cromwell.database.slick.tables

import java.sql.Timestamp
import javax.sql.rowset.serial.SerialClob

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

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255)) // TODO: rename column via liquibase

    def callFullyQualifiedName = column[Option[String]]("CALL_FQN", O.Length(255)) // TODO: rename column via liquibase

    def jobIndex = column[Option[Int]]("JOB_SCATTER_INDEX") // TODO: rename column via liquibase

    def jobAttempt = column[Option[Int]]("JOB_RETRY_ATTEMPT") // TODO: rename column via liquibase

    def metadataKey = column[String]("METADATA_KEY", O.Length(255))

    def metadataValue = column[Option[SerialClob]]("METADATA_VALUE")

    def metadataValueType = column[Option[String]]("METADATA_VALUE_TYPE", O.Length(10))

    def metadataTimestamp = column[Timestamp]("METADATA_TIMESTAMP")

    override def * = (workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt, metadataKey, metadataValue,
      metadataValueType, metadataTimestamp, metadataEntryId.?) <> (MetadataEntry.tupled, MetadataEntry.unapply)

    // TODO: rename index via liquibase
    def ixMetadataEntryWeu = index("METADATA_WORKFLOW_IDX", workflowExecutionUuid, unique = false)
  }

  val metadataEntries = TableQuery[MetadataEntries]

  val metadataEntryIdsAutoInc = metadataEntries returning metadataEntries.map(_.metadataEntryId)

  val metadataEntriesExists = Compiled(metadataEntries.take(1).exists)

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

  def metadataEntryExistsForWorkflowExecutionUuid(workflowId: Rep[String], key: Rep[String]): Rep[Boolean] = {
    metadataEntries.filter( metadataEntry =>
      metadataEntry.workflowExecutionUuid === workflowId &&
      metadataEntry.metadataKey === key &&
      metadataEntry.metadataValue.isDefined
    ).exists
  }

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

  val metadataEntriesForIdRange = Compiled(
    (minMetadataEntryId: Rep[Long], maxMetadataEntryId: Rep[Long],
     startMetadataKey: Rep[String], endMetadataKey: Rep[String],
     nameMetadataKey: Rep[String], statusMetadataKey: Rep[String], submissionMetadataKey: Rep[String],
     parentWorkflowIdMetadataKey: Rep[String], rootWorkflowIdMetadataKey: Rep[String],
     likeLabelMetadataKey: Rep[String]) => {
      for {
        metadataEntry <- metadataEntries
        if metadataEntry.metadataEntryId >= minMetadataEntryId
        if metadataEntry.metadataEntryId <= maxMetadataEntryId
        if metadataEntry.metadataKey === startMetadataKey ||
          metadataEntry.metadataKey === endMetadataKey ||
          metadataEntry.metadataKey === nameMetadataKey ||
          metadataEntry.metadataKey === statusMetadataKey ||
          metadataEntry.metadataKey === submissionMetadataKey ||
          metadataEntry.metadataKey === parentWorkflowIdMetadataKey ||
          metadataEntry.metadataKey === rootWorkflowIdMetadataKey ||
          // http://slick.lightbend.com/doc/3.2.3/queries.html#expressions
          metadataEntry.metadataKey.like(likeLabelMetadataKey ++ ("%": Rep[String]))
        if metadataEntry.callFullyQualifiedName.isEmpty
        if metadataEntry.jobIndex.isEmpty
        if metadataEntry.jobAttempt.isEmpty
      } yield metadataEntry
    }
  )

  val existsMetadataEntriesGreaterThanMetadataEntryId = Compiled(
    (metadataEntryId: Rep[Long]) => {
      val query = for {
        metadataEntry <- metadataEntries
        if metadataEntry.metadataEntryId > metadataEntryId
      } yield metadataEntry.metadataEntryId
      query.exists
    }
  )

  /**
    * Returns metadata entries that are "like" metadataKeys for the specified workflow.
    * If requireEmptyJobKey is true, only workflow level keys are returned, otherwise both workflow and call level
    * keys are returned.
    */
  def metadataEntriesWithKeyConstraints(workflowExecutionUuid: String,
                                        metadataKeysToFilterFor: List[String],
                                        metadataKeysToFilterOut: List[String],
                                        requireEmptyJobKey: Boolean) = {
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeysToFilterFor, metadataKeysToFilterOut)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  }

  /**
    * Returns metadata entries that are "like" metadataKeys for the specified call.
    * If jobAttempt has no value, all metadata keys for all attempts are returned.
    */
  def metadataEntriesForJobWithKeyConstraints(workflowExecutionUuid: String,
                                              metadataKeysToFilterFor: List[String],
                                              metadataKeysToFilterOut: List[String],
                                              callFqn: String,
                                              jobIndex: Option[Int],
                                              jobAttempt: Option[Int]) = {
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeysToFilterFor, metadataKeysToFilterOut)
      if metadataEntry.callFullyQualifiedName === callFqn
      if hasSameIndex(metadataEntry, jobIndex)
      // Assume that every metadata entry for a call should have a non null attempt value
      // Because of that, if the jobAttempt parameter is Some(_), make sure it matches, otherwise take all entries
      // regardless of the attempt
      if (metadataEntry.jobAttempt === jobAttempt) || jobAttempt.isEmpty
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  }

  private[this] def metadataEntryHasMetadataKeysLike(metadataEntry: MetadataEntries,
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
