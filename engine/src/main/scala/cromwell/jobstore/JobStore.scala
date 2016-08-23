package cromwell.jobstore

import cromwell.core.WorkflowId
import wdl4s.TaskOutput

import scala.concurrent.{ExecutionContext, Future}

trait JobStore {
  def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit]
  def readJobResult(jobStoreKey: JobStoreKey, taskOutputs: Seq[TaskOutput])(implicit ec: ExecutionContext): Future[Option[JobResult]]
}
