package cromwell.services.keyvalue

import akka.actor.ActorRef
import cats.data.{NonEmptyList, NonEmptyVector}
import cromwell.core.actor.BatchActor.CommandAndReplyTo
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.EnhancedBatchActor
import cromwell.services.keyvalue.KeyValueServiceActor.{KvFailure, KvPut, KvPutSuccess}

import scala.concurrent.Future
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success}

abstract class KeyValueWriteActor(override val serviceRegistryActor: ActorRef,
                                  flushRate: FiniteDuration,
                                  batchSize: Int
) extends EnhancedBatchActor[CommandAndReplyTo[KvPut]](flushRate, batchSize) {

  override protected def process(data: NonEmptyVector[CommandAndReplyTo[KvPut]]) = instrumentedProcess {
    val processed = processPut(data.map(_.command).toVector)
    processed onComplete {
      case Success(_) =>
        data.toVector.foreach { case CommandAndReplyTo(command, replyTo) =>
          replyTo ! KvPutSuccess(command)
        }
      case Failure(f) =>
        data.toVector.foreach { case CommandAndReplyTo(command, replyTo) =>
          replyTo ! KvFailure(command, f)
        }
    }
    // This method should return how many "operations" have been performed to enable instrumentation of throughput
    // Here we've processed all the KvPuts in "data"
    processed.map(_ => data.length)
  }

  def processPut(put: Vector[KvPut]): Future[Unit]

  // EnhancedBatchActor overrides
  override def receive = enhancedReceive.orElse(super.receive)
  override protected def weightFunction(command: CommandAndReplyTo[KvPut]) = 1
  override protected def instrumentationPath =
    KeyValueServiceActor.InstrumentationPath.concatNel(NonEmptyList.one("write"))
  override protected def instrumentationPrefix = InstrumentationPrefixes.ServicesPrefix
  override def commandToData(snd: ActorRef) = { case put: KvPut =>
    CommandAndReplyTo(put, snd)
  }
}
