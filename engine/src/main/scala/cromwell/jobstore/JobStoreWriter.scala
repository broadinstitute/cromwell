package cromwell.jobstore

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.WorkflowId

import scala.util.{Failure, Success}

/**
  * Singleton actor to coordinate writing job statuses to the database.
  *
  * State: Represents an actor either doing nothing, or currently writing to the database
  * Data: If currently writing, the actor stores pending updates in the data. When one write completes, any further writes are written
  */
case class JobStoreWriter(jsd: JobStoreDatabase) extends LoggingFSM[JobStoreWriterState, JobStoreWriterData] {

  implicit val ec = context.dispatcher

  startWith(Pending, JobStoreWriterData.empty)

  when(Pending) {
    case Event(command: RegisterCompletionMessage, stateData) =>
      val newData = writeNextOperationToDatabase(stateData.withNewOperation(sender, command))
      goto(WritingToDatabase) using newData
  }

  when(WritingToDatabase) {
    case Event(command: RegisterCompletionMessage, stateData) =>
      stay using stateData.withNewOperation(sender, command)
    case Event(WriteComplete, stateData) =>
      val newData = writeNextOperationToDatabase(stateData)
      goto(if (newData.isEmpty) Pending else WritingToDatabase) using newData
  }

  whenUnhandled {
    case Event(someMessage, stateData) =>
      log.error(s"JobStoreWriter: Unexpected message received in state $stateName: $someMessage")
      stay()
  }

  onTransition {
    case (oldState, newState) =>
      log.debug(s"Transitioning from $oldState to $newState")
  }

  def writeNextOperationToDatabase(data: JobStoreWriterData): JobStoreWriterData = {

    val newData = data.rolledOver

    val workflowCompletions = newData.currentOperation collect {
      case (_, RegisterWorkflowCompleted(wfid)) =>  wfid
    }

    val jobCompletions = newData.currentOperation collect {
      case (_, RegisterJobCompleted(jobStoreKey, jobResult)) if !workflowCompletions.contains(jobStoreKey.workflowId) => (jobStoreKey, jobResult)
    }

    if (!(workflowCompletions.isEmpty && jobCompletions.isEmpty)) {
      jsd.writeToDatabase(jobCompletions.toMap, workflowCompletions) onComplete {
        case Success(_) =>
          newData.currentOperation foreach { case (actor, message) => actor ! JobStoreWriteSuccess(message) }
          self ! WriteComplete
        case Failure(reason) =>
          newData.currentOperation foreach { case (actor, message) => actor ! JobStoreWriteFailure(message, reason) }
          self ! WriteComplete
      }
    }

    newData
  }
}

object JobStoreWriter {
  def props(jobStoreDatabase: JobStoreDatabase): Props = Props(new JobStoreWriter(jobStoreDatabase))
}

object JobStoreWriterData {
  def empty = JobStoreWriterData(List.empty, List.empty)
}

case class JobStoreWriterData(currentOperation: List[(ActorRef, RegisterCompletionMessage)], nextOperation: List[(ActorRef, RegisterCompletionMessage)]){
  def isEmpty = nextOperation.isEmpty && currentOperation.isEmpty
  def withNewOperation(sender: ActorRef, command: RegisterCompletionMessage) = this.copy(nextOperation = this.nextOperation :+ (sender, command))
  def rolledOver = JobStoreWriterData(this.nextOperation, List.empty)
}

sealed trait JobStoreWriterState
case object Pending extends JobStoreWriterState
case object WritingToDatabase extends JobStoreWriterState

sealed trait JobStoreWriterMessage
sealed trait RegisterCompletionMessage extends JobStoreWriterMessage
case class RegisterJobCompleted(jobKey: JobStoreKey, jobResult: JobResult) extends RegisterCompletionMessage
case class RegisterWorkflowCompleted(workflowId: WorkflowId) extends RegisterCompletionMessage
case object WriteComplete extends JobStoreWriterMessage

sealed trait JobStoreWriterResponse
case class JobStoreWriteSuccess(writeCommand: RegisterCompletionMessage) extends JobStoreWriterResponse
case class JobStoreWriteFailure(writeCommand: RegisterCompletionMessage, reason: Throwable) extends JobStoreWriterResponse