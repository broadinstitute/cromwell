package cromwell.database.sql

import java.sql.Timestamp

import cromwell.database.sql.tables.WorkflowStoreEntry

import scala.concurrent.{ExecutionContext, Future}

trait WorkflowStoreSqlDatabase {
  this: SqlDatabase =>
  /*
  The following section relates to:
____    __    ____  ______   .______       __  ___  _______  __        ______   ____    __    ____
\   \  /  \  /   / /  __  \  |   _  \     |  |/  / |   ____||  |      /  __  \  \   \  /  \  /   /
 \   \/    \/   / |  |  |  | |  |_)  |    |  '  /  |  |__   |  |     |  |  |  |  \   \/    \/   /
  \            /  |  |  |  | |      /     |    <   |   __|  |  |     |  |  |  |   \            /
   \    /\    /   |  `--'  | |  |\  \----.|  .  \  |  |     |  `----.|  `--'  |    \    /\    /
    \__/  \__/     \______/  | _| `._____||__|\__\ |__|     |_______| \______/      \__/  \__/

     _______.___________.  ______   .______       _______
    /       |           | /  __  \  |   _  \     |   ____|
   |   (----`---|  |----`|  |  |  | |  |_)  |    |  |__
    \   \       |  |     |  |  |  | |      /     |   __|
.----)   |      |  |     |  `--'  | |  |\  \----.|  |____
|_______/       |__|      \______/  | _| `._____||_______|

   */

  /**
    * Set all running workflows to aborting state.
    */
  def setStateToState(fromWorkflowState: String, toWorkflowState: String)
                     (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Set the workflow Id from one state to another.
    *
    * @param workflowExecutionUuid  Id to update or delete.
    * @param workflowStateToDelete1 Delete rows with this state.
    * @param workflowStateToDelete2 Also delete rows with this state.
    * @param workflowStateForUpdate Update other rows to this state.
    * @param ec                     The execution context where to run asynchronous operations.
    * @return None if no rows were matched, Some(true) if a row was deleted or Some(false) if a row was updated.
    */
  def deleteOrUpdateWorkflowToState(workflowExecutionUuid: String,
                                    workflowStateToDelete1: String,
                                    workflowStateToDelete2: String,
                                    workflowStateForUpdate: String)
                                   (implicit ec: ExecutionContext): Future[Option[Boolean]]

  /**
    * Adds the requested WorkflowSourceFiles to the store.
    */
  def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                             (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Retrieves a limited number of workflows which have not already been pulled into the engine and updates their
    * state. NOTE: Rows are returned with the original query state, NOT the updated state.
    *
    * @param limit                      Maximum number of rows to retrieve for this cromwell.
    * @param cromwellId                 The cromwell id identifying  this cromwell.
    * @param heartbeatTimestampTimedOut Any heartbeat earlier than this timestamp is considered free.
    * @param heartbeatTimestampTo       The new timestamp value to write.
    * @param workflowStateFrom          The workflow state that should be searched for, usually 'Submitted'.
    * @param workflowStateTo            The workflow state to change to, usually 'Running'.
    * @param workflowStateExcluded      Workflow states that should not be queried for processing, usually 'On Hold'.
    * @param ec                         The execution context where to run asynchronous operations.
    * @return A list of workflow entries in their original state, NOT the updated state.
    */
  def fetchWorkflowsInState(limit: Int,
                            cromwellId: String,
                            heartbeatTimestampTimedOut: Timestamp,
                            heartbeatTimestampTo: Timestamp,
                            workflowStateFrom: String,
                            workflowStateTo: String,
                            workflowStateExcluded: String)
                           (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]]

  def writeWorkflowHeartbeats(workflowExecutionUuids: Seq[String],
                              heartbeatTimestamp: Timestamp)
                             (implicit ec: ExecutionContext): Future[Int]

  /**
    * Clears out cromwellId and heartbeatTimestamp for all workflow store entries currently assigned
    * the specified cromwellId.
    */
  def releaseWorkflowStoreEntries(cromwellId: String)(implicit ec: ExecutionContext): Future[Int]

  /**
    * Deletes a workflow from the database, returning the number of rows affected.
    */
  def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int]

  def workflowStateCounts(implicit ec: ExecutionContext): Future[Map[String, Int]]

  /**
    * Returns the number of rows updated from one state to another.
    */
  def updateWorkflowState(workflowExecutionUuid: String, fromWorkflowState: String, toWorkflowState: String)
                         (implicit ec: ExecutionContext): Future[Int]

  def findWorkflowsWithAbortRequested(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[String]]

  def findWorkflows(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[String]]

}
