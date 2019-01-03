package cromwell.database.sql

import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState

import scala.concurrent.duration.FiniteDuration
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
  def setAllRunningToAborting()
                             (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Set restarted flags for all running and aborting workflows to true.
    */
  def markRunningAndAbortingAsRestarted(cromwellId: Option[String])
                                       (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Set the workflow to aborting state.
    * @return Some(restarted) if the workflow exists in the store, where restarted is the value of its restarted flag
    *         None if the workflow does not exist in the store
    */
  def setToAborting(workflowId: String)
                   (implicit ec: ExecutionContext): Future[Option[Boolean]]

  /**
    * Adds the requested WorkflowSourceFiles to the store.
    */
  def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                             (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Retrieves up to limit workflows which have not already been pulled into the engine and updates their state.
    * NOTE: Rows are returned with the query state, NOT the update state.
    */
  def fetchStartableWorkflows(limit: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                             (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]]

  def writeWorkflowHeartbeats(workflowExecutionUuids: Set[String])(implicit ec: ExecutionContext): Future[Int]

  /**
    * Clears out cromwellId and heartbeatTimestamp for all workflow store entries currently assigned
    * the specified cromwellId.
    */
  def releaseWorkflowStoreEntries(cromwellId: String)(implicit ec: ExecutionContext): Future[Unit]

  /**
    * Deletes a workflow from the database, returning the number of rows affected.
    */
  def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int]
  
  def stats(implicit ec: ExecutionContext): Future[Map[WorkflowStoreState, Int]]

  def setOnHoldToSubmitted(workflowId: String)(implicit ec: ExecutionContext): Future[Unit]
}
