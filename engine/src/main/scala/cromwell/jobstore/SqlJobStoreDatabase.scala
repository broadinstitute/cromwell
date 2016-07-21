package cromwell.jobstore

import cromwell.core.WorkflowId
import cromwell.database.CromwellDatabase

import scala.concurrent.{Future, ExecutionContext}

class SqlJobStoreDatabase extends JobStoreDatabase with CromwellDatabase {
  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = ???

  override def readJobResult(jobStoreKey: JobStoreKey)(implicit ec: ExecutionContext): Future[Option[JobResult]] = ???
}
