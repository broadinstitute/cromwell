package cromwell.jobstore

import akka.actor.{LoggingFSM, Props}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.actor.BatchingDbWriter._
import cromwell.core.actor.{BatchingDbWriter, BatchingDbWriterActor}
import cromwell.jobstore.JobStore.{JobCompletion, WorkflowCompletion}
import cromwell.jobstore.JobStoreActor._

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}


case class JobStoreWriterActor(jsd: JobStore, dbBatchSize: Int, override val dbFlushRate: FiniteDuration) extends LoggingFSM[BatchingDbWriterState, BatchingDbWriter.BatchingDbWriterData] with BatchingDbWriterActor {

  implicit val ec = context.dispatcher

  log.info(s"JobStoreWriterActor configured to write to the database with batch size {} and flush rate {}.", dbBatchSize, dbFlushRate)

  startWith(WaitingToWrite, NoData)

  when(WaitingToWrite) {
    case Event(command: JobStoreWriterCommand, curData) =>
      curData.addData(CommandAndReplyTo(command, sender)) match {
        case newData: HasData[_] if newData.length >= dbBatchSize => goto(WritingToDb) using newData
        case newData => stay() using newData
      }
    case Event(ScheduledFlushToDb, curData) =>
      log.debug("Initiating periodic job store flush to DB")
      goto(WritingToDb) using curData
  }

  when(WritingToDb) {
    case Event(ScheduledFlushToDb, _) => stay
    case Event(command: JobStoreWriterCommand, curData) => stay using curData.addData(CommandAndReplyTo(command, sender))
    case Event(FlushBatchToDb, NoData) =>
      log.debug("Attempted job store flush to DB but had nothing to write")
      goto(WaitingToWrite)
    case Event(FlushBatchToDb, HasData(data)) =>
      log.debug("Flushing {} job store commands to the DB", data.length)
      val completions = data.toVector.collect({ case CommandAndReplyTo(c: JobStoreWriterCommand, _) => c.completion })

      if (completions.nonEmpty) {
        val workflowCompletions = completions collect { case w: WorkflowCompletion => w }
        val completedWorkflowIds = workflowCompletions map { _.workflowId } toSet
        // Filter job completions that also have a corresponding workflow completion; these would just be
        // immediately deleted anyway.
        val jobCompletions = completions.toList collect { case j: JobCompletion if !completedWorkflowIds.contains(j.key.workflowId) => j }

        jsd.writeToDatabase(workflowCompletions, jobCompletions, dbBatchSize) onComplete {
          case Success(_) =>
            data map { case CommandAndReplyTo(c: JobStoreWriterCommand, r) => r ! JobStoreWriteSuccess(c) }
            self ! DbWriteComplete
          case Failure(regerts) =>
            log.error("Failed to properly job store entries to database", regerts)
            data map { case CommandAndReplyTo(_, r) => r ! JobStoreWriteFailure(regerts) }
            self ! DbWriteComplete
        }
      }
      stay using NoData
    case Event(DbWriteComplete, _) =>
      log.debug("Flush of job store commands complete")
      goto(WaitingToWrite)
  }
}

object JobStoreWriterActor {

  def props(jobStoreDatabase: JobStore, dbBatchSize: Int, dbFlushRate: FiniteDuration): Props = Props(new JobStoreWriterActor(jobStoreDatabase, dbBatchSize, dbFlushRate)).withDispatcher(EngineDispatcher)
}
