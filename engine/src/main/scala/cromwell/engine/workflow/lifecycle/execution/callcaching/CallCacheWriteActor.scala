package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, Props}
import cats.data.NonEmptyVector
import cats.instances.list._
import cats.instances.tuple._
import cats.syntax.foldable._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.actor.BatchActor
import cromwell.core.actor.BatchActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache.CallCacheHashBundle
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheWriteActor.SaveCallCacheHashes

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps

case class CallCacheWriteActor(callCache: CallCache) extends BatchActor[CommandAndReplyTo[SaveCallCacheHashes]](CallCacheWriteActor.dbFlushRate, CallCacheWriteActor.dbBatchSize) {

  def commandToData(snd: ActorRef): PartialFunction[Any, CommandAndReplyTo[SaveCallCacheHashes]] = {
    case command: SaveCallCacheHashes => CommandAndReplyTo(command, snd)
  }

  override protected def weightFunction(command: CommandAndReplyTo[SaveCallCacheHashes]) = 1

  override protected def process(data: NonEmptyVector[CommandAndReplyTo[SaveCallCacheHashes]]) = {
    log.debug("Flushing {} call cache hashes sets to the DB", data.length)

    //     Collect all the bundles of hashes that should be written and all the senders which should be informed of
    //     success or failure.
    val (bundles, replyTos) = data.toList.foldMap { case CommandAndReplyTo(s: SaveCallCacheHashes, r: ActorRef) => (List(s.bundle), List(r)) }
    if (bundles.nonEmpty) {
      val futureMessage = callCache.addToCache(bundles, batchSize) map { _ => CallCacheWriteSuccess } recover { case t => CallCacheWriteFailure(t) }
      futureMessage map { message =>
        replyTos foreach { _ ! message }
      }
      futureMessage.map(_ => data.length)
    } else Future.successful(0)
  }
}

object CallCacheWriteActor {
  def props(callCache: CallCache): Props = Props(CallCacheWriteActor(callCache)).withDispatcher(EngineDispatcher)

  case class SaveCallCacheHashes(bundle: CallCacheHashBundle)

  val dbBatchSize = 100
  val dbFlushRate = 3 seconds
}

sealed trait CallCacheWriteResponse
case object CallCacheWriteSuccess extends CallCacheWriteResponse
case class CallCacheWriteFailure(t: Throwable) extends CallCacheWriteResponse
