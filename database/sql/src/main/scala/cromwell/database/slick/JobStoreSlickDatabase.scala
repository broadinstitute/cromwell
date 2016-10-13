package cromwell.database.slick

import cats.instances.future._
import cats.syntax.functor._
import cromwell.database.sql.JobStoreSqlDatabase
import cromwell.database.sql.joins.JobStoreJoin

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

trait JobStoreSlickDatabase extends JobStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  private def addJobStore(jobStoreJoin: JobStoreJoin)(implicit ec: ExecutionContext): DBIO[Unit] = {
    for {
      jobStoreEntryId <- dataAccess.jobStoreEntryIdsAutoInc += jobStoreJoin.jobStoreEntry
      _ <- dataAccess.jobStoreSimpletonEntryIdsAutoInc ++=
        jobStoreJoin.jobStoreSimpletonEntries.map(_.copy(jobStoreEntryId = Option(jobStoreEntryId)))
    } yield ()
  }

  override def addJobStores(jobStoreJoins: Seq[JobStoreJoin])
                           (implicit ec: ExecutionContext): Future[Unit] = {
    val action = DBIO.sequence(jobStoreJoins map addJobStore)
    runTransaction(action) void
  }

  override def queryJobStores(workflowExecutionUuid: String, callFqn: String, jobScatterIndex: Int,
                              jobScatterAttempt: Int)(implicit ec: ExecutionContext):
  Future[Option[JobStoreJoin]] = {

    val action = for {
      jobStoreEntryOption <- dataAccess.
        jobStoreEntriesForJobKey((workflowExecutionUuid, callFqn, jobScatterIndex, jobScatterAttempt)).result.headOption
      jobStoreSimpletonEntries <- jobStoreEntryOption match {
        case Some(jobStoreEntry) =>
          dataAccess.jobStoreSimpletonEntriesForJobStoreEntryId(jobStoreEntry.jobStoreEntryId.get).result
        case _ => DBIO.successful(Seq.empty)
      }
    } yield jobStoreEntryOption.map(JobStoreJoin(_, jobStoreSimpletonEntries))

    runTransaction(action)
  }

  override def removeJobStores(workflowExecutionUuids: Seq[String])
                              (implicit ec: ExecutionContext): Future[Seq[Int]] = {
    val actions = workflowExecutionUuids map {
      workflowExecutionUuid => dataAccess.jobStoreEntriesForWorkflowExecutionUuid(workflowExecutionUuid).delete
    }
    runTransaction(DBIO.sequence(actions))
  }
}
