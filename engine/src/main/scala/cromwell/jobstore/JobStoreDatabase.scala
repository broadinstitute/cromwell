package cromwell.jobstore

import cromwell.core.WorkflowId

import scala.concurrent.{ExecutionContext, Future}

trait JobStoreDatabase {
  def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit]
}
