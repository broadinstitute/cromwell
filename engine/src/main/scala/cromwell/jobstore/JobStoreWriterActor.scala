package cromwell.jobstore

import akka.actor.{ActorRef, Props}
import cats.data.{NonEmptyList, NonEmptyVector}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.LoadConfig
import cromwell.core.actor.BatchActor._
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.engine.workflow.workflowstore.WorkflowStoreAccess
import cromwell.jobstore.JobStore.{JobCompletion, WorkflowCompletion}
import cromwell.jobstore.JobStoreActor._
import cromwell.services.EnhancedBatchActor

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}


case class JobStoreWriterActor(jsd: JobStore,
                               override val batchSize: Int,
                               override val flushRate: FiniteDuration,
                               serviceRegistryActor: ActorRef,
                               threshold: Int,
                               workflowStoreAccess: WorkflowStoreAccess
                              )
  extends EnhancedBatchActor[CommandAndReplyTo[JobStoreWriterCommand]](flushRate, batchSize) {

  override protected def process(nonEmptyData: NonEmptyVector[CommandAndReplyTo[JobStoreWriterCommand]]): Future[Int] = instrumentedProcess {
    val data = nonEmptyData.toVector
    log.debug("Flushing {} job store commands to the DB", data.length)
    val completions = data.collect({ case CommandAndReplyTo(c: JobStoreWriterCommand, _) => c.completion })

    if (completions.nonEmpty) {
      val workflowCompletions = completions collect { case w: WorkflowCompletion => w }
      val completedWorkflowIds = workflowCompletions map { _.workflowId } toSet
      // Filter job completions that also have a corresponding workflow completion; these would just be
      // immediately deleted anyway.
      val jobCompletions = completions.toList collect { case j: JobCompletion if !completedWorkflowIds.contains(j.key.workflowId) => j }
      val jobStoreAction: Future[Unit] = jsd.writeToDatabase(workflowCompletions, jobCompletions, batchSize)
      val workflowStoreAction: Future[List[Int]] = Future.sequence {
        completedWorkflowIds.map(workflowStoreAccess.deleteFromStore(_)).toList
      }

      val combinedAction: Future[Unit] = for {
        _ <- jobStoreAction
        _ <- workflowStoreAction
      } yield ()

      combinedAction onComplete {
        case Success(_) =>
          data foreach { case CommandAndReplyTo(c: JobStoreWriterCommand, r) => r ! JobStoreWriteSuccess(c) }
        case Failure(regerts) =>
          log.error(regerts, "Failed to write job store entries to database")
          data foreach { case CommandAndReplyTo(_, r) => r ! JobStoreWriteFailure(regerts) }
      }

      combinedAction.map(_ => 1)
    } else Future.successful(0)
  }

  // EnhancedBatchActor overrides
  override def receive = enhancedReceive.orElse(super.receive)
  override protected def weightFunction(command: CommandAndReplyTo[JobStoreWriterCommand]) = 1
  override protected def instrumentationPath = NonEmptyList.of("store", "write")
  override protected def instrumentationPrefix = InstrumentationPrefixes.JobPrefix
  override def commandToData(snd: ActorRef): PartialFunction[Any, CommandAndReplyTo[JobStoreWriterCommand]] = {
    case command: JobStoreWriterCommand => CommandAndReplyTo(command, snd)
  }
}

object JobStoreWriterActor {
  def props(jobStoreDatabase: JobStore,
            dbBatchSize: Int,
            dbFlushRate: FiniteDuration,
            registryActor: ActorRef,
            workflowStoreAccess: WorkflowStoreAccess): Props = {
    Props(new JobStoreWriterActor(jobStoreDatabase, dbBatchSize, dbFlushRate, registryActor, LoadConfig.JobStoreWriteThreshold, workflowStoreAccess)).withDispatcher(EngineDispatcher)
  }
}
