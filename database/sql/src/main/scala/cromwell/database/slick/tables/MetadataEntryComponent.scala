package cromwell.database.slick.tables

import java.sql.Timestamp

import javax.sql.rowset.serial.SerialClob
import cromwell.database.sql.tables.{InformationSchemaEntry, MetadataEntry}
import slick.jdbc.GetResult

trait MetadataEntryComponent {

  this: DriverComponent with WorkflowMetadataSummaryEntryComponent =>

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

    def workflowExecutionUuid =
      column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255)) // TODO: rename column via liquibase

    def callFullyQualifiedName = column[Option[String]]("CALL_FQN", O.Length(255)) // TODO: rename column via liquibase

    def jobIndex = column[Option[Int]]("JOB_SCATTER_INDEX") // TODO: rename column via liquibase

    def jobAttempt = column[Option[Int]]("JOB_RETRY_ATTEMPT") // TODO: rename column via liquibase

    def metadataKey = column[String]("METADATA_KEY", O.Length(255))

    def metadataValue = column[Option[SerialClob]]("METADATA_VALUE")

    def metadataValueType = column[Option[String]]("METADATA_VALUE_TYPE", O.Length(10))

    def metadataTimestamp = column[Timestamp]("METADATA_TIMESTAMP")

    override def * = (workflowExecutionUuid,
                      callFullyQualifiedName,
                      jobIndex,
                      jobAttempt,
                      metadataKey,
                      metadataValue,
                      metadataValueType,
                      metadataTimestamp,
                      metadataEntryId.?
    ) <> (MetadataEntry.tupled, MetadataEntry.unapply)

    // TODO: rename index via liquibase
    def ixMetadataEntryWeu = index("METADATA_WORKFLOW_IDX", workflowExecutionUuid, unique = false)
  }

  val metadataEntries = TableQuery[MetadataEntries]

  val metadataEntryIdsAutoInc = metadataEntries returning metadataEntries.map(_.metadataEntryId)

  val metadataEntriesExists = Compiled(metadataEntries.take(1).exists)

  val metadataEntriesForWorkflowExecutionUuid = Compiled((workflowExecutionUuid: Rep[String]) =>
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val metadataEntriesForWorkflowSortedById = Compiled((workflowExecutionUuid: Rep[String]) =>
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).sortBy(_.metadataEntryId)
  )

  val countMetadataEntriesForWorkflowExecutionUuid =
    Compiled((rootWorkflowId: Rep[String], expandSubWorkflows: Rep[Boolean]) =>
      {
        val targetWorkflowIds = for {
          summary <- workflowMetadataSummaryEntries
          // Uses `IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU`, `UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU`
          if summary.workflowExecutionUuid === rootWorkflowId || ((summary.rootWorkflowExecutionUuid === rootWorkflowId) && expandSubWorkflows)
        } yield summary.workflowExecutionUuid

        for {
          metadata <- metadataEntries
          if metadata.workflowExecutionUuid in targetWorkflowIds // Uses `METADATA_WORKFLOW_IDX`
        } yield metadata
      }.size
    )

  val metadataEntryExistsForWorkflowExecutionUuid = Compiled((workflowExecutionUuid: Rep[String]) =>
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield metadataEntry).exists
  )

  def metadataEntryExistsForWorkflowExecutionUuid(workflowId: Rep[String], key: Rep[String]): Rep[Boolean] =
    metadataEntries
      .filter(metadataEntry =>
        metadataEntry.workflowExecutionUuid === workflowId &&
          metadataEntry.metadataKey === key &&
          metadataEntry.metadataValue.isDefined
      )
      .exists

  val metadataEntriesForWorkflowExecutionUuidAndMetadataKey =
    Compiled((workflowExecutionUuid: Rep[String], metadataKey: Rep[String]) =>
      (for {
        metadataEntry <- metadataEntries
        if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
        if metadataEntry.metadataKey === metadataKey
        if metadataEntry.callFullyQualifiedName.isEmpty
        if metadataEntry.jobIndex.isEmpty
        if metadataEntry.jobAttempt.isEmpty
      } yield metadataEntry).sortBy(_.metadataTimestamp)
    )

  val countMetadataEntriesForWorkflowExecutionUuidAndMetadataKey =
    Compiled((rootWorkflowId: Rep[String], metadataKey: Rep[String], expandSubWorkflows: Rep[Boolean]) =>
      {
        val targetWorkflowIds = for {
          summary <- workflowMetadataSummaryEntries
          // Uses `IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU`, `UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU`
          if summary.workflowExecutionUuid === rootWorkflowId || ((summary.rootWorkflowExecutionUuid === rootWorkflowId) && expandSubWorkflows)
        } yield summary.workflowExecutionUuid

        for {
          metadata <- metadataEntries
          if metadata.workflowExecutionUuid in targetWorkflowIds // Uses `METADATA_WORKFLOW_IDX`
          if metadata.metadataKey === metadataKey
          if metadata.callFullyQualifiedName.isEmpty
          if metadata.jobIndex.isEmpty
          if metadata.jobAttempt.isEmpty
        } yield metadata
      }.size
    )

  val metadataEntriesForJobKey = Compiled(
    (workflowExecutionUuid: Rep[String],
     callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Option[Int]]
    ) =>
      (for {
        metadataEntry <- metadataEntries
        if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
        if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
        if hasSameIndex(metadataEntry, jobIndex)
        if hasSameAttempt(metadataEntry, jobAttempt)
      } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val countMetadataEntriesForJobKey = Compiled(
    (rootWorkflowId: Rep[String],
     callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Option[Int]],
     expandSubWorkflows: Rep[Boolean]
    ) =>
      {
        val targetWorkflowIds = for {
          summary <- workflowMetadataSummaryEntries
          // Uses `IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU`, `UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU`
          if summary.workflowExecutionUuid === rootWorkflowId || ((summary.rootWorkflowExecutionUuid === rootWorkflowId) && expandSubWorkflows)
        } yield summary.workflowExecutionUuid

        for {
          metadata <- metadataEntries
          if metadata.workflowExecutionUuid in targetWorkflowIds // Uses `METADATA_WORKFLOW_IDX`
          if metadata.callFullyQualifiedName === callFullyQualifiedName
          if hasSameIndex(metadata, jobIndex)
          if hasSameAttempt(metadata, jobAttempt)
        } yield metadata
      }.size
  )

  val metadataEntriesForJobKeyAndMetadataKey = Compiled(
    (workflowExecutionUuid: Rep[String],
     metadataKey: Rep[String],
     callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Option[Int]]
    ) =>
      (for {
        metadataEntry <- metadataEntries
        if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
        if metadataEntry.metadataKey === metadataKey
        if metadataEntry.callFullyQualifiedName === callFullyQualifiedName
        if hasSameIndex(metadataEntry, jobIndex)
        if hasSameAttempt(metadataEntry, jobAttempt)
      } yield metadataEntry).sortBy(_.metadataTimestamp)
  )

  val countMetadataEntriesForJobKeyAndMetadataKey = Compiled(
    (rootWorkflowId: Rep[String],
     metadataKey: Rep[String],
     callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Option[Int]],
     jobAttempt: Rep[Option[Int]],
     expandSubWorkflows: Rep[Boolean]
    ) =>
      {
        val targetWorkflowIds = for {
          summary <- workflowMetadataSummaryEntries
          // Uses `IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU`, `UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU`
          if summary.workflowExecutionUuid === rootWorkflowId || ((summary.rootWorkflowExecutionUuid === rootWorkflowId) && expandSubWorkflows)
        } yield summary.workflowExecutionUuid

        for {
          metadata <- metadataEntries
          if metadata.workflowExecutionUuid in targetWorkflowIds // Uses `METADATA_WORKFLOW_IDX`
          if metadata.metadataKey === metadataKey
          if metadata.callFullyQualifiedName === callFullyQualifiedName
          if hasSameIndex(metadata, jobIndex)
          if hasSameAttempt(metadata, jobAttempt)
        } yield metadata
      }.size
  )

  val metadataEntriesForIdRange = Compiled { (minMetadataEntryId: Rep[Long], maxMetadataEntryId: Rep[Long]) =>
    for {
      metadataEntry <- metadataEntries
      if metadataEntry.metadataEntryId >= minMetadataEntryId
      if metadataEntry.metadataEntryId <= maxMetadataEntryId
    } yield metadataEntry
  }

  /**
    * Returns metadata entries that are "like" metadataKeys for the specified workflow.
    * If requireEmptyJobKey is true, only workflow level keys are returned, otherwise both workflow and call level
    * keys are returned.
    */
  def metadataEntriesWithKeyConstraints(workflowExecutionUuid: String,
                                        metadataKeysToFilterFor: List[String],
                                        metadataKeysToFilterOut: List[String],
                                        requireEmptyJobKey: Boolean
  ) =
    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid === workflowExecutionUuid
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeysToFilterFor, metadataKeysToFilterOut)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry).sortBy(_.metadataTimestamp)

  /**
    * Counts metadata entries that are "like" metadataKeys for the specified workflow.
    * If requireEmptyJobKey is true, only workflow level keys are counted, otherwise both workflow and call level
    * keys are counted.
    */
  def countMetadataEntriesWithKeyConstraints(rootWorkflowId: String,
                                             metadataKeysToFilterFor: List[String],
                                             metadataKeysToFilterOut: List[String],
                                             requireEmptyJobKey: Boolean,
                                             expandSubWorkflows: Boolean
  ) = {

    val targetWorkflowIds = for {
      summary <- workflowMetadataSummaryEntries
      // Uses `IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU`, `UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU`
      if summary.workflowExecutionUuid === rootWorkflowId || ((summary.rootWorkflowExecutionUuid === rootWorkflowId) && expandSubWorkflows)
    } yield summary.workflowExecutionUuid

    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid in targetWorkflowIds
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeysToFilterFor, metadataKeysToFilterOut)
      if metadataEntryHasEmptyJobKey(metadataEntry, requireEmptyJobKey)
    } yield metadataEntry).size
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
                                              jobAttempt: Option[Int]
  ) =
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

  /**
    * Counts metadata entries that are "like" metadataKeys for the specified call.
    * If jobAttempt has no value, all metadata keys for all attempts are counted.
    */
  def countMetadataEntriesForJobWithKeyConstraints(rootWorkflowId: String,
                                                   metadataKeysToFilterFor: List[String],
                                                   metadataKeysToFilterOut: List[String],
                                                   callFqn: String,
                                                   jobIndex: Option[Int],
                                                   jobAttempt: Option[Int],
                                                   expandSubWorkflows: Boolean
  ) = {

    val targetWorkflowIds = for {
      summary <- workflowMetadataSummaryEntries
      // Uses `IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU`, `UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU`
      if summary.workflowExecutionUuid === rootWorkflowId || ((summary.rootWorkflowExecutionUuid === rootWorkflowId) && expandSubWorkflows)
    } yield summary.workflowExecutionUuid

    (for {
      metadataEntry <- metadataEntries
      if metadataEntry.workflowExecutionUuid in targetWorkflowIds
      if metadataEntryHasMetadataKeysLike(metadataEntry, metadataKeysToFilterFor, metadataKeysToFilterOut)
      if metadataEntry.callFullyQualifiedName === callFqn
      if hasSameIndex(metadataEntry, jobIndex)
      // Assume that every metadata entry for a call should have a non null attempt value
      // Because of that, if the jobAttempt parameter is Some(_), make sure it matches, otherwise take all entries
      // regardless of the attempt
      if (metadataEntry.jobAttempt === jobAttempt) || jobAttempt.isEmpty
    } yield metadataEntry).size
  }

  def metadataTableSizeInformation() = {
    val query =
      sql"""
           |SELECT DATA_LENGTH, INDEX_LENGTH, DATA_FREE
           |FROM information_schema.tables
           |WHERE TABLE_NAME = 'METADATA_ENTRY'
         """.stripMargin
    query
      .as[InformationSchemaEntry](rconv = GetResult { r =>
        InformationSchemaEntry(r.<<, r.<<, r.<<)
      })
      .headOption
  }

  def failedJobsMetadataWithWorkflowId(rootWorkflowId: String, isPostgres: Boolean) = {
    val getMetadataEntryResult = GetResult { r =>
      MetadataEntry(r.<<,
                    r.<<,
                    r.<<,
                    r.<<,
                    r.<<,
                    r.nextClobOption().map(clob => new SerialClob(clob)),
                    r.<<,
                    r.<<,
                    r.<<
      )
    }

    def dbIdentifierWrapper(identifier: String, isPostgres: Boolean) =
      if (isPostgres) s"${'"'}$identifier${'"'}" else identifier

    def evaluateMetadataValue(isPostgres: Boolean, colName: String): String =
      if (isPostgres) s"convert_from(lo_get(${colName}::oid), 'UTF8')" else colName

    def attemptAndIndexSelectStatement(callFqn: String,
                                       scatterIndex: String,
                                       retryAttempt: String,
                                       variablePrefix: String
    ): String =
      s"SELECT ${callFqn}, MAX(COALESCE(${scatterIndex}, 0)) as ${variablePrefix}Scatter, MAX(COALESCE(${retryAttempt}, 0)) AS ${variablePrefix}Retry"

    val workflowUuid = dbIdentifierWrapper("WORKFLOW_EXECUTION_UUID", isPostgres)
    val callFqn = dbIdentifierWrapper("CALL_FQN", isPostgres)
    val scatterIndex = dbIdentifierWrapper("JOB_SCATTER_INDEX", isPostgres)
    val retryAttempt = dbIdentifierWrapper("JOB_RETRY_ATTEMPT", isPostgres)
    val metadataKey = dbIdentifierWrapper("METADATA_KEY", isPostgres)
    val metadataValueType = dbIdentifierWrapper("METADATA_VALUE_TYPE", isPostgres)
    val metadataTimestamp = dbIdentifierWrapper("METADATA_TIMESTAMP", isPostgres)
    val metadataJournalId = dbIdentifierWrapper("METADATA_JOURNAL_ID", isPostgres)
    val rootUuid = dbIdentifierWrapper("ROOT_WORKFLOW_EXECUTION_UUID", isPostgres)
    val metadataValue = dbIdentifierWrapper("METADATA_VALUE", isPostgres)
    val metadataEntry = dbIdentifierWrapper("METADATA_ENTRY", isPostgres)
    val wmse = dbIdentifierWrapper("WORKFLOW_METADATA_SUMMARY_ENTRY", isPostgres)
    val resultSetColumnNames =
      s"me.${workflowUuid}, me.${callFqn}, me.${scatterIndex}, me.${retryAttempt}, me.${metadataKey}, me.${metadataValue}, me.${metadataValueType}, me.${metadataTimestamp}, me.${metadataJournalId}"

    val query =
      sql"""
      SELECT #${resultSetColumnNames}
      FROM #${metadataEntry} me
      INNER JOIN (
       #${attemptAndIndexSelectStatement(callFqn, scatterIndex, retryAttempt, "failed")}
        FROM #${metadataEntry} me
        INNER JOIN #${wmse} wmse
        ON wmse.#${workflowUuid} = me.#${workflowUuid}
        WHERE (wmse.#${rootUuid} = $rootWorkflowId OR wmse.#${workflowUuid} = $rootWorkflowId)
        AND (me.#${metadataKey} in ('executionStatus', 'backendStatus') AND #${evaluateMetadataValue(isPostgres,
                                                                                                     metadataValue
        )} = 'Failed')
        GROUP BY #${callFqn}, #${metadataValue} 
        HAVING #${evaluateMetadataValue(isPostgres, metadataValue)} = 'Failed'
      ) AS failedCalls
      ON me.#${callFqn} = failedCalls.#${callFqn}
      INNER JOIN (
        #${attemptAndIndexSelectStatement(callFqn, scatterIndex, retryAttempt, "max")}
        FROM #${metadataEntry} me
        INNER JOIN #${wmse} wmse
        ON wmse.#${workflowUuid} = me.#${workflowUuid}
        WHERE (wmse.#${rootUuid} = $rootWorkflowId OR wmse.#${workflowUuid} = $rootWorkflowId)
        AND #${callFqn} IS NOT NULL
        GROUP BY #${callFqn}
      ) maxCalls
      ON me.#${callFqn} = maxCalls.#${callFqn}
      LEFT JOIN (
        SELECT DISTINCT #${callFqn}
        FROM #${metadataEntry} me
        INNER JOIN #${wmse} wmse
        ON wmse.#${workflowUuid} = me.#${workflowUuid}
        WHERE (wmse.#${rootUuid} = $rootWorkflowId OR wmse.#${workflowUuid} = $rootWorkflowId)
        AND me.#${metadataKey} = 'subWorkflowId'
        GROUP BY #${callFqn}
      ) AS avoidedCalls
      ON me.#${callFqn} = avoidedCalls.#${callFqn}
      INNER JOIN #${wmse} wmse
      ON wmse.#${workflowUuid} = me.#${workflowUuid}
      WHERE avoidedCalls.#${callFqn} IS NULL
      AND COALESCE(me.#${scatterIndex}, 0) = maxCalls.maxScatter
      AND COALESCE(me.#${retryAttempt}, 0) = maxCalls.maxRetry
      AND failedCalls.failedScatter = maxCalls.maxScatter
      AND failedCalls.failedRetry = maxCalls.maxRetry
      GROUP BY #${resultSetColumnNames}
      HAVING me.#${workflowUuid} IN (
        SELECT DISTINCT wmse.#${workflowUuid}
        FROM #${wmse} wmse
        WHERE wmse.#${rootUuid} = $rootWorkflowId OR wmse.#${workflowUuid} = $rootWorkflowId
      )
    """

    query.as(getMetadataEntryResult)
  }

  private[this] def metadataEntryHasMetadataKeysLike(metadataEntry: MetadataEntries,
                                                     metadataKeysToFilterFor: List[String],
                                                     metadataKeysToFilterOut: List[String]
  ): Rep[Boolean] = {

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

  private[this] def hasSameIndex(metadataEntry: MetadataEntries, jobIndex: Rep[Option[Int]]) =
    (metadataEntry.jobIndex.isEmpty && jobIndex.isEmpty) || (metadataEntry.jobIndex === jobIndex)

  private[this] def hasSameAttempt(metadataEntry: MetadataEntries, jobAttempt: Rep[Option[Int]]) =
    (metadataEntry.jobAttempt.isEmpty && jobAttempt.isEmpty) || (metadataEntry.jobAttempt === jobAttempt)

  private[this] def metadataEntryHasEmptyJobKey(metadataEntry: MetadataEntries,
                                                requireEmptyJobKey: Rep[Boolean]
  ): Rep[Boolean] =
    !requireEmptyJobKey ||
      (metadataEntry.callFullyQualifiedName.isEmpty &&
        metadataEntry.jobIndex.isEmpty &&
        metadataEntry.jobAttempt.isEmpty)
}
