package cromwell.database.sql

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

  def updateWorkflowState(queryWorkflowState: String, updateWorkflowState: String)
                         (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Adds the requested WorkflowSourceFiles to the store.
    */
  def addWorkflowStoreEntries(workflowStoreEntries: Iterable[WorkflowStoreEntry])
                             (implicit ec: ExecutionContext): Future[Unit]

  /**
    * Retrieves up to limit workflows which have not already been pulled into the engine and updates their state.
    * NOTE: Rows are returned with the query state, NOT the update state.
    */
  def queryWorkflowStoreEntries(limit: Int, queryWorkflowState: String, updateWorkflowState: String)
                               (implicit ec: ExecutionContext): Future[Seq[WorkflowStoreEntry]]

  /**
    * Deletes a workflow from the database, returning the number of rows affected.
    */
  def removeWorkflowStoreEntry(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Int]
}
