package cromwell.jobstore

import akka.actor.{Actor, Props}
import cromwell.core.WorkflowId
import cromwell.jobstore.JobStoreActor.{JobStoreReaderCommand, JobStoreWriterCommand}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage


/**
  * Joins the service registry API to the JobStoreReaderActor and JobStoreWriterActor.
  *
  * This level of indirection is a tiny bit awkward but allows the database to be injected.
  */
class JobStoreActor extends Actor {
  // TODO: Replace with a real database, probably from config.
  val database = FilesystemJobStoreDatabase
  val jobStoreWriterActor = context.actorOf(JobStoreWriterActor.props(database))
  val jobStoreReaderActor = context.actorOf(JobStoreReaderActor.props(database))

  override def receive: Receive = {
    case command: JobStoreWriterCommand => jobStoreWriterActor.tell(command, sender())
    case command: JobStoreReaderCommand => jobStoreReaderActor.tell(command, sender())
  }
}

object JobStoreActor {
  sealed trait JobStoreCommand

  sealed trait JobStoreWriterCommand extends JobStoreCommand
  case class RegisterJobCompleted(jobKey: JobStoreKey, jobResult: JobResult) extends JobStoreWriterCommand
  case class RegisterWorkflowCompleted(workflowId: WorkflowId) extends JobStoreWriterCommand

  sealed trait JobStoreWriterResponse
  case class JobStoreWriteSuccess(originalCommand: JobStoreWriterCommand) extends JobStoreWriterResponse
  case class JobStoreWriteFailure(originalCommand: JobStoreWriterCommand, reason: Throwable) extends JobStoreWriterResponse

  sealed trait JobStoreReaderCommand extends JobStoreCommand
  /**
    * Message to query the JobStoreReaderActor, asks whether the specified job has already been completed.
    */
  case class QueryJobCompletion(jobKey: JobStoreKey) extends JobStoreReaderCommand

  sealed trait JobStoreReaderResponse
  /**
    * Message which indicates that a job has already completed, and contains the results of the job
    */
  case class JobComplete(jobResult: JobResult) extends JobStoreReaderResponse
  /**
    * Indicates that the job has not been completed yet. Makes no statement about whether the job is
    * running versus unstarted or (maybe?) doesn't even exist!
    */
  case object JobNotComplete extends JobStoreReaderResponse

  case class JobStoreReadFailure(reason: Throwable) extends JobStoreReaderResponse

  def props = Props(new JobStoreActor)
}
