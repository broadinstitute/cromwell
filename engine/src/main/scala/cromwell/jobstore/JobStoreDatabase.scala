package cromwell.jobstore

import cromwell.core.WorkflowId

import scala.concurrent.Future

trait JobStoreDatabase {
  def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId]): Future[Unit]
}
