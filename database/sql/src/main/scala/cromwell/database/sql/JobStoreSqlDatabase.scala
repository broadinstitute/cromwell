package cromwell.database.sql

import cromwell.database.sql.joins.JobStoreJoin

import scala.concurrent.{ExecutionContext, Future}

trait JobStoreSqlDatabase {
  this: SqlDatabase =>

  def addJobStores(jobStoreJoins: Seq[JobStoreJoin])(implicit ec: ExecutionContext): Future[Unit]

  def queryJobStores(workflowExecutionUuid: String, callFqn: String, jobScatterIndex: Int, jobScatterAttempt: Int)
                    (implicit ec: ExecutionContext): Future[Option[JobStoreJoin]]

  def removeJobStores(workflowExecutionUuids: Seq[String])(implicit ec: ExecutionContext): Future[Seq[Int]]
}
