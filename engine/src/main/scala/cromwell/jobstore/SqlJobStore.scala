package cromwell.jobstore

import cromwell.core.ExecutionIndex._
import cromwell.core.WorkflowId
import cromwell.core.simpleton.WdlValueSimpleton._
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import cromwell.database.sql.JobStoreSqlDatabase
import cromwell.database.sql.tables.{JobStoreEntry, JobStoreResultSimpletonEntry}
import wdl4s.TaskOutput
import wdl4s.types._
import wdl4s.values.WdlPrimitive

import scala.concurrent.{ExecutionContext, Future}

class SqlJobStore(sqlDatabase: JobStoreSqlDatabase) extends JobStore {
  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = {
    for {
      _ <- sqlDatabase.addJobStoreEntries(jobCompletions map toDatabase)
      _ <- sqlDatabase.removeJobResultsForWorkflows(workflowCompletions map {_.toString})
    } yield ()
  }

  private def toDatabase(jobCompletion: (JobStoreKey, JobResult)): (JobStoreEntry, Int => Iterable[JobStoreResultSimpletonEntry]) = {
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
        val simpletons: Int => Iterable[JobStoreResultSimpletonEntry] = fk => jobOutputs.mapValues(_.wdlValue).simplify.map {
          s => JobStoreResultSimpletonEntry(s.simpletonKey, s.simpletonValue.valueString, s.simpletonValue.wdlType.toWdlString, fk) }
        (entry, simpletons)
      case (key, JobResultFailure(returnCode, throwable, retryable)) =>
        val entry = JobStoreEntry(
          key.workflowId.toString,
          key.callFqn,
          key.index.fromIndex,
          key.attempt,
          jobSuccessful = false,
          returnCode,
          Some(throwable.getMessage),
          Some(retryable))
        val simpletons: Int => Iterable[JobStoreResultSimpletonEntry] = _ => List.empty
        (entry, simpletons)
    }
  }

  private def toSimpleton(entry: JobStoreResultSimpletonEntry): WdlValueSimpleton = {
    val wdlType: WdlType = entry.wdlType match {
      case "String" => WdlStringType
      case "Int" => WdlIntegerType
      case "Float" => WdlFloatType
      case "Boolean" => WdlBooleanType
      case _ => throw new RuntimeException(s"$entry: unrecognized WDL type: ${entry.wdlType}")
    }
    WdlValueSimpleton(entry.simpletonKey, wdlType.coerceRawValue(entry.simpletonValue).get.asInstanceOf[WdlPrimitive])
  }

  override def readJobResult(jobStoreKey: JobStoreKey, taskOutputs: Seq[TaskOutput])(implicit ec: ExecutionContext): Future[Option[JobResult]] = {
    sqlDatabase.fetchJobResult(jobStoreKey.workflowId.toString, jobStoreKey.callFqn, jobStoreKey.index.fromIndex, jobStoreKey.attempt) map {
      _ map { case (entry, simpletonEntries) =>
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
