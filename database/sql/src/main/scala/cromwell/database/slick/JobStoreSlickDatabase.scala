package cromwell.database.slick

import cromwell.database.sql.JobStoreSqlDatabase
import cromwell.database.sql.joins.JobStoreJoin
import cromwell.database.sql.tables.JobStoreSimpletonEntry

import scala.concurrent.{ExecutionContext, Future}

trait JobStoreSlickDatabase extends JobStoreSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

  override def addJobStores(jobStoreJoins: Seq[JobStoreJoin], batchSize: Int)
                           (implicit ec: ExecutionContext): Future[Unit] = {

    def assignJobStoreIdsToSimpletons(jobStoreIds: Seq[Int]): Seq[JobStoreSimpletonEntry] = {
      val simpletonsByJobStoreEntry = jobStoreJoins map { _.jobStoreSimpletonEntries }
      val jobStoreIdsAndSimpletons = jobStoreIds.toList zip simpletonsByJobStoreEntry
      jobStoreIdsAndSimpletons flatMap { case (id, simpletons) => simpletons.map(_.copy(jobStoreEntryId = Option(id))) }
    }

    val jobStoreEntries = jobStoreJoins map { _.jobStoreEntry }

    val action = for {
      jobStoreIds <- dataAccess.jobStoreEntryIdsAutoInc ++= jobStoreEntries.toList
      simpletons = assignJobStoreIdsToSimpletons(jobStoreIds)
      simpletonInserts = simpletons.grouped(batchSize) map { dataAccess.jobStoreSimpletonEntries ++= _ }
      _ <- DBIO.sequence(simpletonInserts)
    } yield ()

    runTransaction(action)
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
