package cromwell.jobstore

import cromwell.Simpletons._
import cromwell.core.ExecutionIndex._
import cromwell.core.WorkflowId
import cromwell.core.simpleton.WdlValueBuilder
import cromwell.core.simpleton.WdlValueSimpleton._
import cromwell.database.sql.JobStoreSqlDatabase
import cromwell.database.sql.joins.JobStoreJoin
import cromwell.database.sql.tables.{JobStoreEntry, JobStoreSimpletonEntry}
import wdl4s.TaskOutput

import scala.concurrent.{ExecutionContext, Future}

class SqlJobStore(sqlDatabase: JobStoreSqlDatabase) extends JobStore {
  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = {
    for {
      _ <- sqlDatabase.addJobStores(jobCompletions.toSeq map toDatabase)
      _ <- sqlDatabase.removeJobStores(workflowCompletions.map(_.toString))
    } yield ()
  }

  private def toDatabase(jobCompletion: (JobStoreKey, JobResult)): JobStoreJoin = {
    jobCompletion match {
      case (key, JobResultSuccess(returnCode, jobOutputs)) =>
        val entry = JobStoreEntry(
          key.workflowId.toString,
          key.callFqn,
          key.index.fromIndex,
          key.attempt,
          jobSuccessful = true,
          returnCode,
          None,
          None)
        val jobStoreResultSimpletons =
          jobOutputs.mapValues(_.wdlValue).simplify.map {
            wdlValueSimpleton => JobStoreSimpletonEntry(
              wdlValueSimpleton.simpletonKey, wdlValueSimpleton.simpletonValue.valueString,
              wdlValueSimpleton.simpletonValue.wdlType.toWdlString)
          }
        JobStoreJoin(entry, jobStoreResultSimpletons.toSeq)
      case (key, JobResultFailure(returnCode, throwable, retryable)) =>
        val entry = JobStoreEntry(
          key.workflowId.toString,
          key.callFqn,
          key.index.fromIndex,
          key.attempt,
          jobSuccessful = false,
          returnCode,
          Option(throwable.getMessage),
          Option(retryable))
        JobStoreJoin(entry, Seq.empty)
    }
  }

  override def readJobResult(jobStoreKey: JobStoreKey, taskOutputs: Seq[TaskOutput])(implicit ec: ExecutionContext): Future[Option[JobResult]] = {
    sqlDatabase.queryJobStores(jobStoreKey.workflowId.toString, jobStoreKey.callFqn, jobStoreKey.index.fromIndex,
      jobStoreKey.attempt) map {
      _ map { case JobStoreJoin(entry, simpletonEntries) =>
        entry match {
          case JobStoreEntry(_, _, _, _, true, returnCode, None, None, _) =>
            val simpletons = simpletonEntries map toSimpleton
            val jobOutputs = WdlValueBuilder.toJobOutputs(taskOutputs, simpletons)
            JobResultSuccess(returnCode, jobOutputs)
          case JobStoreEntry(_, _, _, _, false, returnCode, Some(exceptionMessage), Some(retryable), _) =>
            JobResultFailure(returnCode, new Exception(exceptionMessage), retryable)
          case bad =>
            throw new Exception(s"Invalid contents of JobStore table: $bad")
        }
      }
    }
  }
}
