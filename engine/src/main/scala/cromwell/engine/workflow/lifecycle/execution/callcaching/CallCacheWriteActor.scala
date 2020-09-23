package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, Props}
import cats.data.{NonEmptyList, NonEmptyVector}
import cats.implicits._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.LoadConfig
import cromwell.core.actor.BatchActor._
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache.CallCacheHashBundle
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheWriteActor.SaveCallCacheHashes
import cromwell.services.EnhancedBatchActor

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps

case class CallCacheWriteActor(callCache: CallCache, serviceRegistryActor: ActorRef, threshold: Int)
  extends EnhancedBatchActor[CommandAndReplyTo[SaveCallCacheHashes]](
    CallCacheWriteActor.dbFlushRate,
    CallCacheWriteActor.dbBatchSize) {

  override protected def process(data: NonEmptyVector[CommandAndReplyTo[SaveCallCacheHashes]]) = instrumentedProcess {
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

  // EnhancedBatchActor overrides
  override def receive = enhancedReceive.orElse(super.receive)
  override protected def weightFunction(command: CommandAndReplyTo[SaveCallCacheHashes]) = 1
  override protected def instrumentationPath = NonEmptyList.of("callcaching", "write")
  override protected def instrumentationPrefix = InstrumentationPrefixes.JobPrefix
  def commandToData(snd: ActorRef): PartialFunction[Any, CommandAndReplyTo[SaveCallCacheHashes]] = {
    case command: SaveCallCacheHashes => CommandAndReplyTo(command, snd)
  }
}

object CallCacheWriteActor {
  def props(callCache: CallCache, registryActor: ActorRef): Props = {
    Props(CallCacheWriteActor(callCache, registryActor, LoadConfig.CallCacheWriteThreshold)).withDispatcher(EngineDispatcher)
  }

  case class SaveCallCacheHashes(bundle: CallCacheHashBundle)

  val dbBatchSize = 100
  val dbFlushRate = 3 seconds
}

sealed trait CallCacheWriteResponse
case object CallCacheWriteSuccess extends CallCacheWriteResponse
case class CallCacheWriteFailure(t: Throwable) extends CallCacheWriteResponse
