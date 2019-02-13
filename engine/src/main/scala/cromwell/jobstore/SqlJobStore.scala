package cromwell.jobstore

import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromwell.Simpletons._
import cromwell.backend.async.JobAlreadyFailedInJobStore
import cromwell.core.ExecutionIndex._
import cromwell.core.simpleton.WomValueBuilder
import cromwell.core.simpleton.WomValueSimpleton._
import cromwell.database.sql.EngineSqlDatabase
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.JobStoreJoin
import cromwell.database.sql.tables.{JobStoreEntry, JobStoreSimpletonEntry}
import cromwell.jobstore.JobStore.{JobCompletion, WorkflowCompletion}
import org.slf4j.LoggerFactory
import wom.graph.GraphNodePort.OutputPort

import scala.concurrent.{ExecutionContext, Future}

class SqlJobStore(sqlDatabase: EngineSqlDatabase) extends JobStore {
  val log = LoggerFactory.getLogger(classOf[SqlJobStore])

  override def writeToDatabase(workflowCompletions: Seq[WorkflowCompletion], jobCompletions: Seq[JobCompletion], batchSize: Int)(implicit ec: ExecutionContext): Future[Unit] = {
    val completedWorkflowIds = workflowCompletions.toList.map(_.workflowId.toString)
    for {
      _ <- sqlDatabase.addJobStores(jobCompletions map toDatabase, batchSize)
      _ <- completedWorkflowIds traverse sqlDatabase.removeWorkflowStoreEntry
      _ <- completedWorkflowIds traverse sqlDatabase.removeDockerHashStoreEntries
      _ <- sqlDatabase.removeJobStores(completedWorkflowIds)
    } yield ()
  }

  private def toDatabase(jobCompletion: JobCompletion): JobStoreJoin = {
    jobCompletion match {
      case JobCompletion(key, JobResultSuccess(returnCode, jobOutputs)) =>
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
          jobOutputs.outputs.simplify.map {
            womValueSimpleton => JobStoreSimpletonEntry(
              womValueSimpleton.simpletonKey, womValueSimpleton.simpletonValue.valueString.toClobOption,
              womValueSimpleton.simpletonValue.womType.stableName)
          }
        JobStoreJoin(entry, jobStoreResultSimpletons.toSeq)
      case JobCompletion(key, JobResultFailure(returnCode, throwable, retryable)) =>
        val entry = JobStoreEntry(
          key.workflowId.toString,
          key.callFqn,
          key.index.fromIndex,
          key.attempt,
          jobSuccessful = false,
          returnCode,
          Option(throwable.getMessage).toClobOption,
          Option(retryable))
        JobStoreJoin(entry, Seq.empty)
    }
  }

  override def readJobResult(jobStoreKey: JobStoreKey, taskOutputs: Seq[OutputPort])(implicit ec: ExecutionContext): Future[Option[JobResult]] = {
    sqlDatabase.queryJobStores(jobStoreKey.workflowId.toString, jobStoreKey.callFqn, jobStoreKey.index.fromIndex,
      jobStoreKey.attempt) map {
      _ map { case JobStoreJoin(entry, simpletonEntries) =>
        entry match {
          case JobStoreEntry(_, _, _, _, true, returnCode, None, None, _) =>
            val simpletons = simpletonEntries map toSimpleton
            val jobOutputs = WomValueBuilder.toJobOutputs(taskOutputs, simpletons)
            JobResultSuccess(returnCode, jobOutputs)
          case JobStoreEntry(_, _, _, _, false, returnCode, Some(_), Some(retryable), _) =>
            JobResultFailure(returnCode,
              JobAlreadyFailedInJobStore(jobStoreKey.tag, entry.exceptionMessage.toRawString),
              retryable)
          case bad =>
            throw new Exception(s"Invalid contents of JobStore table: $bad")
        }
      }
    }
  }
}
