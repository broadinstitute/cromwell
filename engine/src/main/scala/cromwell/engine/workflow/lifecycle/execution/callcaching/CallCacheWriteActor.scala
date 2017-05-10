package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, LoggingFSM, Props}
import cats.instances.list._
import cats.instances.tuple._
import cats.syntax.foldable._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.core.actor.BatchingDbWriter
import cromwell.core.actor.BatchingDbWriter._
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.SucceededResponseData
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache.CallCacheHashBundle
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheWriteActor.{SaveCallCacheHashes, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CallCacheHashes

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

case class CallCacheWriteActor(callCache: CallCache) extends LoggingFSM[BatchingDbWriterState, BatchingDbWriter.BatchingDbWriterData] {

  implicit val ec: ExecutionContext = context.dispatcher

  startWith(WaitingToWrite, NoData)

  override def preStart(): Unit = {
    context.system.scheduler.schedule(0.seconds, dbFlushRate, self, ScheduledFlushToDb)
    super.preStart()
  }

  when(WaitingToWrite) {
    case Event(command: SaveCallCacheHashes, curData) =>
      curData.addData(CommandAndReplyTo(command, sender)) match {
        case newData: HasData[_] if newData.length >= dbBatchSize => goto(WritingToDb) using newData
        case newData => stay() using newData
      }
    case Event(ScheduledFlushToDb, curData) =>
      log.debug("Initiating periodic call cache flush to DB")
      goto(WritingToDb)
  }

  when(WritingToDb) {
    case Event(ScheduledFlushToDb, curData) => stay
    case Event(command: SaveCallCacheHashes, curData) => stay using curData.addData(CommandAndReplyTo(command, sender))
    case Event(FlushBatchToDb, NoData) =>
      log.debug("Attempted call cache hash set flush to DB but had nothing to write")
      goto(WaitingToWrite)
    case Event(FlushBatchToDb, HasData(data)) =>
      log.debug("Flushing {} call cache hashes sets to the DB", data.length)

      // Collect all the bundles of hashes that should be written and all the senders which should be informed of
      // success or failure.
      val (bundles, replyTos) = data.foldMap { case CommandAndReplyTo(s: SaveCallCacheHashes, r: ActorRef) => (List(s.bundle), List(r)) }
      if (bundles.nonEmpty) {
        val futureMessage = callCache.addToCache(bundles, dbBatchSize) map { _ => CallCacheWriteSuccess } recover { case t => CallCacheWriteFailure(t) }
        futureMessage map { message =>
          replyTos foreach { _ ! message }
          self ! DbWriteComplete
        }
      }
      stay using NoData
    case Event(DbWriteComplete, curData) =>
      log.debug("Flush of cache data complete")
      goto(WaitingToWrite)
  }

  onTransition {
    case WaitingToWrite -> WritingToDb => self ! FlushBatchToDb
  }
}

object CallCacheWriteActor {
  def props(callCache: CallCache): Props = Props(CallCacheWriteActor(callCache)).withDispatcher(EngineDispatcher)

  case class SaveCallCacheHashes(workflowId: WorkflowId, hashes: CallCacheHashes, data: SucceededResponseData) {
    def bundle = CallCacheHashBundle(workflowId, hashes, data.successResponse)
  }

  val dbBatchSize = 100
  val dbFlushRate = 3 seconds
}

sealed trait CallCacheWriteResponse
case object CallCacheWriteSuccess extends CallCacheWriteResponse
case class CallCacheWriteFailure(t: Throwable) extends CallCacheWriteResponse
