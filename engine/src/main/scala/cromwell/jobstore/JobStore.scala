package cromwell.jobstore

import cromwell.core.WorkflowId
import cromwell.jobstore.JobStore.{JobCompletion, WorkflowCompletion}
import wom.graph.GraphNodePort.OutputPort

import scala.concurrent.{ExecutionContext, Future}

trait JobStore {
  def writeToDatabase(workflowCompletions: Seq[WorkflowCompletion], jobCompletions: Seq[JobCompletion], batchSize: Int)(implicit ec: ExecutionContext): Future[Unit]
  def readJobResult(jobStoreKey: JobStoreKey, taskOutputs: Seq[OutputPort])(implicit ec: ExecutionContext): Future[Option[JobResult]]
}

object JobStore {
  sealed trait Completion
  case class WorkflowCompletion(workflowId: WorkflowId) extends Completion
  case class JobCompletion(key: JobStoreKey, result: JobResult) extends Completion
}
