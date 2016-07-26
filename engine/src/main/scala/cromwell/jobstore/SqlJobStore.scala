package cromwell.jobstore

import cromwell.core.{JobOutputs, WorkflowId}
import cromwell.database.sql.JobStoreSqlDatabase
import cromwell.database.sql.tables.JobStoreEntry

import scala.concurrent.{ExecutionContext, Future}
import cromwell.core.ExecutionIndex._
import spray.json._

class SqlJobStore(sqlDatabase: JobStoreSqlDatabase) extends JobStore {
  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = {
    for {
      _ <- sqlDatabase.addJobStoreEntries(jobCompletions map toJobStoreEntry)
      _ <- sqlDatabase.removeJobResultsForWorkflows(workflowCompletions map {_.toString})
    } yield ()
  }

  private def toJobStoreEntry(jobCompletion: (JobStoreKey, JobResult)): JobStoreEntry = {
    import JobResultJsonFormatter._
    jobCompletion match {
      case (key, JobResultSuccess(returnCode, jobResult)) => JobStoreEntry(
        key.workflowId.toString,
        key.callFqn,
        key.index.fromIndex,
        key.attempt,
        jobSuccessful = true,
        returnCode,
        Some(jobResult.toJson.toString),
        None)
      case (key, JobResultFailure(returnCode, throwable)) => JobStoreEntry(
        key.workflowId.toString,
        key.callFqn,
        key.index.fromIndex,
        key.attempt,
        jobSuccessful = true,
        returnCode,
        None,
        Some(throwable.getMessage))
    }
  }

  override def readJobResult(jobStoreKey: JobStoreKey)(implicit ec: ExecutionContext): Future[Option[JobResult]] = {
    sqlDatabase.fetchJobResult(jobStoreKey.workflowId.toString, jobStoreKey.callFqn, jobStoreKey.index.fromIndex, jobStoreKey.attempt) map { _ map { toJobResult } }
  }

  private def toJobResult(jobStoreEntry: JobStoreEntry) = {
    import JobResultJsonFormatter._

    jobStoreEntry match {
      case JobStoreEntry(_, _, _, _, true, returnCode, Some(jobResult), _, _) =>
        JobResultSuccess(returnCode, jobResult.parseJson.convertTo[JobOutputs])
      case JobStoreEntry(_, _, _, _, false, returnCode, _, Some(exceptionMessage), _) =>
        JobResultFailure(returnCode, new Exception(exceptionMessage))
      case bad => throw new Exception(s"Invalid contents of JobStore table: $bad")
    }
  }
}
