package cromwell.jobstore

import akka.actor.LoggingFSM
import cromwell.core.{JobKey, WorkflowId}

class JobStoreWriter extends LoggingFSM[JobStoreWriterState, JobStoreWriterData] {

}

object JobStoreWriter {

}

case class JobStoreWriterData(jobCompletionsList: Map[JobKey, JobResult], workflowCompletionsList: List[WorkflowId])

sealed trait JobStoreWriterState
case object NothinDoin extends JobStoreWriterState
case object WritingToDatabase extends JobStoreWriterState

sealed trait JobStoreWriterCommand
case class RegisterJobCompleted(jobKey: JobKey, jobResult: JobResult)
case class RegisterWorkflowCompleted(workflowId: WorkflowId)