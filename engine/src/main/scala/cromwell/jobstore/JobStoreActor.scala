package cromwell.jobstore

import akka.actor.{Actor, ActorLogging, Props}
import cats.data.NonEmptyList
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.jobstore.JobStore.{Completion, JobCompletion, WorkflowCompletion}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import wdl4s.wdl.TaskOutput

import scala.concurrent.duration._
import scala.language.postfixOps
/**
  * Joins the service registry API to the JobStoreReaderActor and JobStoreWriterActor.
  *
  * This level of indirection is a tiny bit awkward but allows the database to be injected.
  */
class JobStoreActor(jobStore: JobStore, dbBatchSize: Int, dbFlushRate: FiniteDuration) extends Actor with ActorLogging with GracefulShutdownHelper {
  import JobStoreActor._
  val jobStoreWriterActor = context.actorOf(JobStoreWriterActor.props(jobStore, dbBatchSize, dbFlushRate), "JobStoreWriterActor")
  val jobStoreReaderActor = context.actorOf(JobStoreReaderActor.props(jobStore), "JobStoreReaderActor")

  override def receive: Receive = {
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(jobStoreWriterActor))
    case command: JobStoreWriterCommand => jobStoreWriterActor.tell(command, sender())
    case command: JobStoreReaderCommand => jobStoreReaderActor.tell(command, sender())
  }
}

object JobStoreActor {
  sealed trait JobStoreCommand

  sealed trait JobStoreWriterCommand extends JobStoreCommand {
    def completion: Completion
  }
  case class RegisterJobCompleted(jobKey: JobStoreKey, jobResult: JobResult) extends JobStoreWriterCommand {
    override def completion = JobCompletion(jobKey, jobResult)
  }
  case class RegisterWorkflowCompleted(workflowId: WorkflowId) extends JobStoreWriterCommand {
    override def completion = WorkflowCompletion(workflowId)
  }

  sealed trait JobStoreWriterResponse
  case class JobStoreWriteSuccess(originalCommand: JobStoreWriterCommand) extends JobStoreWriterResponse
  case class JobStoreWriteFailure(reason: Throwable) extends JobStoreWriterResponse

  sealed trait JobStoreReaderCommand extends JobStoreCommand
  /**
    * Message to query the JobStoreReaderActor, asks whether the specified job has already been completed.
    */
  case class QueryJobCompletion(jobKey: JobStoreKey, taskOutputs: Seq[TaskOutput]) extends JobStoreReaderCommand

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

  def props(database: JobStore) = Props(new JobStoreActor(database, dbBatchSize, dbFlushRate)).withDispatcher(EngineDispatcher)

  val dbFlushRate = 1 second

  // `dbBatchSize` applies to simpletons only.  Batching job store entry writes while returning the inserted IDs works at the
  // Slick API level, but didn't actually batch the SQL.  Unfortunately these IDs are required to assign into the simpletons.
  val dbBatchSize = 1000
}
