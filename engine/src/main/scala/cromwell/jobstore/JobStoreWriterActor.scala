package cromwell.jobstore

import akka.actor.{ActorRef, Props}
import cats.data.NonEmptyVector
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.actor.BatchActor._
import cromwell.core.actor.BatchActor
import cromwell.jobstore.JobStore.{JobCompletion, WorkflowCompletion}
import cromwell.jobstore.JobStoreActor._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}


case class JobStoreWriterActor(jsd: JobStore, override val batchSize: Int, override val flushRate: FiniteDuration) extends 
  BatchActor[CommandAndReplyTo[JobStoreWriterCommand]](flushRate, batchSize) {

  override protected def weightFunction(command: CommandAndReplyTo[JobStoreWriterCommand]) = 1

  def commandToData(snd: ActorRef): PartialFunction[Any, CommandAndReplyTo[JobStoreWriterCommand]] = {
    case command: JobStoreWriterCommand => CommandAndReplyTo(command, snd)
  }

  override protected def process(nonEmptyData: NonEmptyVector[CommandAndReplyTo[JobStoreWriterCommand]]) = {
    val data = nonEmptyData.toVector
    log.debug("Flushing {} job store commands to the DB", data.length)
    val completions = data.collect({ case CommandAndReplyTo(c: JobStoreWriterCommand, _) => c.completion })

    if (completions.nonEmpty) {
      val workflowCompletions = completions collect { case w: WorkflowCompletion => w }
      val completedWorkflowIds = workflowCompletions map { _.workflowId } toSet
      // Filter job completions that also have a corresponding workflow completion; these would just be
      // immediately deleted anyway.
      val jobCompletions = completions.toList collect { case j: JobCompletion if !completedWorkflowIds.contains(j.key.workflowId) => j }
      val databaseAction = jsd.writeToDatabase(workflowCompletions, jobCompletions, batchSize)

      databaseAction onComplete {
        case Success(_) =>
          data foreach { case CommandAndReplyTo(c: JobStoreWriterCommand, r) => r ! JobStoreWriteSuccess(c) }
        case Failure(regerts) =>
          log.error("Failed to properly job store entries to database", regerts)
          data foreach { case CommandAndReplyTo(_, r) => r ! JobStoreWriteFailure(regerts) }
      }

      databaseAction.map(_ => 1)
    } else Future.successful(0)
  }
}

object JobStoreWriterActor {

  def props(jobStoreDatabase: JobStore, dbBatchSize: Int, dbFlushRate: FiniteDuration): Props = Props(new JobStoreWriterActor(jobStoreDatabase, dbBatchSize, dbFlushRate)).withDispatcher(EngineDispatcher)
}
