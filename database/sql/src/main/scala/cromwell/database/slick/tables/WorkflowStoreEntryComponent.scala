package cromwell.database.slick.tables

import java.sql.Timestamp
import javax.sql.rowset.serial.{SerialBlob, SerialClob}
import cromwell.database.sql.tables.WorkflowStoreEntry

trait WorkflowStoreEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class WorkflowStoreEntries(tag: Tag) extends Table[WorkflowStoreEntry](tag, "WORKFLOW_STORE_ENTRY") {
    def workflowStoreEntryId = column[Int]("WORKFLOW_STORE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255))

    def workflowRoot = column[Option[String]]("WORKFLOW_ROOT", O.Length(100))

    def workflowType = column[Option[String]]("WORKFLOW_TYPE", O.Length(30))

    def workflowTypeVersion = column[Option[String]]("WORKFLOW_TYPE_VERSION", O.Length(255))

    def workflowDefinition = column[Option[SerialClob]]("WORKFLOW_DEFINITION")

    def workflowUrl = column[Option[String]]("WORKFLOW_URL", O.Length(2000))

    def workflowInputs = column[Option[SerialClob]]("WORKFLOW_INPUTS")

    def workflowOptions = column[Option[SerialClob]]("WORKFLOW_OPTIONS")

    def customLabels = column[SerialClob]("CUSTOM_LABELS")

    def workflowState = column[String]("WORKFLOW_STATE", O.Length(20))

    def submissionTime = column[Timestamp]("SUBMISSION_TIME")

    def importsZip = column[Option[SerialBlob]]("IMPORTS_ZIP")

    def cromwellId = column[Option[String]]("CROMWELL_ID", O.Length(100))

    def heartbeatTimestamp = column[Option[Timestamp]]("HEARTBEAT_TIMESTAMP")

    def hogGroup = column[Option[String]]("HOG_GROUP", O.Length(100))

    override def * = (workflowExecutionUuid, workflowDefinition, workflowUrl, workflowRoot, workflowType, workflowTypeVersion, workflowInputs, workflowOptions, workflowState,
      submissionTime, importsZip, customLabels, cromwellId, heartbeatTimestamp, hogGroup, workflowStoreEntryId.?) <> ((WorkflowStoreEntry.apply _).tupled, WorkflowStoreEntry.unapply)

    def ucWorkflowStoreEntryWeu = index("UC_WORKFLOW_STORE_ENTRY_WEU", workflowExecutionUuid, unique = true)

    def ixWorkflowStoreEntryWs = index("IX_WORKFLOW_STORE_ENTRY_WS", workflowState, unique = false)
  }

  protected val workflowStoreEntries = TableQuery[WorkflowStoreEntries]

  val workflowStoreEntryIdsAutoInc = workflowStoreEntries returning workflowStoreEntries.map(_.workflowStoreEntryId)

  /**
    * Useful for finding the workflow store for a given workflow execution UUID
    */
  val workflowStoreEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowStoreEntry
  )

  val heartbeatForWorkflowStoreEntry = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowStoreEntry.heartbeatTimestamp
  )

  /**
   * Return hog group with the lowest count of actively running workflows, that has nonzero startable workflows
   */
  def getHogGroupWithLowestRunningWfs(heartbeatTimestampTimedOut: Timestamp,
                                      excludeWorkflowState: String,
                                      excludedGroups: Set[String]): Query[Rep[Option[String]], Option[String], Seq] = {
    val startableWorkflows = for {
      row <- workflowStoreEntries
      /*
        This looks for:
        1) Workflows with no heartbeat (newly submitted or from a cleanly shut down Cromwell).
        2) Workflows with old heartbeats, presumably abandoned by a defunct Cromwell.
        3) Workflows not in "OnHold" state
        4) Workflows that don't belong to hog groups in excludedGroups
      */
      if (row.heartbeatTimestamp.isEmpty || row.heartbeatTimestamp < heartbeatTimestampTimedOut) &&
        (row.workflowState =!= excludeWorkflowState) &&
        !(row.hogGroup inSet excludedGroups)
    } yield row

    // calculates the count of startable workflows per hog group and oldest submission timestamp of workflow
    // present in that hog group
    val numOfStartableWfsByHogGroup: Query[
      (Rep[Option[String]], Rep[Int], Rep[Option[Timestamp]]),
      (Option[String], Int, Option[Timestamp]),
      Seq
    ] = startableWorkflows
      .groupBy(_.hogGroup)
      .map { case (hogGroupName, workflows) => (hogGroupName, workflows.length, workflows.map(_.submissionTime).min) }
      .sortBy(_._2.asc)

    val totalWorkflows = for {
      row <- workflowStoreEntries
      /*
        This looks for:
        1) Workflows not in "OnHold" state
        2) Workflows that don't belong to hog groups in excludedGroups
      */
      if row.workflowState =!= excludeWorkflowState &&
        !(row.hogGroup inSet excludedGroups)
    } yield row

    // calculates the count of total workflows per hog group
    val totalWorkflowsByHogGroup: Query[(Rep[Option[String]], Rep[Int]), (Option[String], Int), Seq] = totalWorkflows
      .groupBy(_.hogGroup)
      .map { case (hogGroupName, workflows) => (hogGroupName, workflows.length) }
      .sortBy(_._2.asc)

    // calculates the number of actively running workflows for each hog group. If a hog group
    // has all it's workflows that are either actively running or in "OnHold" status it is not
    // included in the list. Hog groups that have no workflows actively running return count as 0
    val wfsRunningPerHogGroup: Query[
      (Rep[Option[String]], Rep[Int], Rep[Option[Timestamp]]),
      (Option[String], Int, Option[Timestamp]),
      Seq
    ] = for {
      (hog_group, workflows_ct) <- totalWorkflowsByHogGroup
      (startable_hog_group, startable_workflows_ct, oldest_submission_time) <- numOfStartableWfsByHogGroup if hog_group === startable_hog_group
    } yield (hog_group, workflows_ct - startable_workflows_ct, oldest_submission_time)

    // sort the above calculated result set first by the count of actively running workflows, then by hog group with
    // oldest submission timestamp and then sort it alphabetically by hog group name. Then take the first row of
    // the result and return the hog group name.
    wfsRunningPerHogGroup.sortBy {
      case (hogGroupName, running_wf_ct, oldest_submission_time) => (running_wf_ct.asc, oldest_submission_time, hogGroupName)
    }.take(1).map(_._1)
  }

  /**
   * Returns up to "limit" startable workflows, sorted by submission time, that belong to
   * given hog group and are not in "OnHold" status.
   */
  val fetchStartableWfsForHogGroup = Compiled(
    (limit: ConstColumn[Long],
     heartbeatTimestampTimedOut: ConstColumn[Timestamp],
     excludeWorkflowState: Rep[String],
     hogGroup: Rep[Option[String]]) => {

      val workflowsToStart = for {
        row <- workflowStoreEntries
        /*
          This looks for:
          1) Workflows with no heartbeat (newly submitted or from a cleanly shut down Cromwell).
          2) Workflows with old heartbeats, presumably abandoned by a defunct Cromwell.
          3) Workflows not in "OnHold" state
          4) Workflows that belong to included hog group
        */
        if (row.heartbeatTimestamp.isEmpty || row.heartbeatTimestamp < heartbeatTimestampTimedOut) &&
          (row.workflowState =!= excludeWorkflowState) &&
          (row.hogGroup === hogGroup)
      } yield row

      /*
        This is a "query for update", meaning rows are locked such that readers are blocked since we will
        do an update subsequent to this select in the same transaction that we know will impact those readers.
       */
      workflowsToStart.forUpdate.sortBy(_.submissionTime.asc).take(limit)
    }
  )

  /**
    * Useful for counting workflows in a given state.
    */
  val workflowStoreStats = Compiled(
    for {
      (state, entry) <- workflowStoreEntries groupBy (_.workflowState)
    } yield state -> entry.size
  )

  /**
    * Useful for updating the relevant fields of a workflow store entry when a workflow is picked up for processing.
    */
  val workflowStoreFieldsForPickup = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      row <- workflowStoreEntries
      if row.workflowExecutionUuid === workflowExecutionUuid
    } yield (row.workflowState, row.cromwellId, row.heartbeatTimestamp)
  )

  /**
    * Useful for clearing out cromwellId and heartbeatTimestamp on an orderly Cromwell shutdown.
    */
  val releaseWorkflowStoreEntries = Compiled(
    (cromwellId: Rep[String]) => for {
      row <- workflowStoreEntries
      if row.cromwellId === cromwellId
    } yield (row.cromwellId, row.heartbeatTimestamp)
  )

  /**
    * Useful for updating state for all entries matching a given state
    */
  val workflowStateForWorkflowState = Compiled(
    (workflowState: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowState === workflowState
    } yield workflowStoreEntry.workflowState
  )

  /**
    * Useful for updating a given workflow to a new state
    */
  val workflowStateForWorkflowExecutionUUid = Compiled(
    (workflowId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowId
    } yield workflowStoreEntry.workflowState
  )

  /**
    * Useful for updating a given workflow to a 'Submitted' state when it's currently 'On Hold'
    */
  val workflowStateForWorkflowExecutionUUidAndWorkflowState = Compiled(
    (workflowId: Rep[String], workflowState: Rep[String]) => {
      for {
        workflowStoreEntry <- workflowStoreEntries
        if workflowStoreEntry.workflowExecutionUuid === workflowId
        if workflowStoreEntry.workflowState === workflowState
      } yield workflowStoreEntry.workflowState
    }
  )

  /**
    * Useful for deleting a given workflow to a 'Submitted' state when it's currently 'On Hold' or 'Submitted'
    */
  val workflowStoreEntryForWorkflowExecutionUUidAndWorkflowStates = Compiled(
    (workflowId: Rep[String],
     workflowStateOr1: Rep[String],
     workflowStateOr2: Rep[String]
    ) => {
      for {
        workflowStoreEntry <- workflowStoreEntries
        if workflowStoreEntry.workflowExecutionUuid === workflowId
        if workflowStoreEntry.workflowState === workflowStateOr1 ||
          workflowStoreEntry.workflowState === workflowStateOr2
      } yield workflowStoreEntry
    }
  )

  // Find workflows running on a given Cromwell instance with abort requested:
  val findWorkflowsWithAbortRequested = Compiled(
    (cromwellId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowState === "Aborting" && workflowStoreEntry.cromwellId === cromwellId
    } yield workflowStoreEntry.workflowExecutionUuid
  )

  // Find workflows running on a given Cromwell instance:
  val findWorkflows = Compiled(
    (cromwellId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.cromwellId === cromwellId
    } yield workflowStoreEntry.workflowExecutionUuid
  )

  val checkExists = Compiled(
    (workflowId: Rep[String]) => (for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowId
    } yield 1)
  )
}
