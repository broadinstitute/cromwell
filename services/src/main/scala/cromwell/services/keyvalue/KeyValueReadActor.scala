package cromwell.services.keyvalue

import akka.actor.ActorRef
import cats.data.NonEmptyList
import cromwell.core.actor.BatchActor.CommandAndReplyTo
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.EnhancedThrottlerActor
import cromwell.services.keyvalue.KeyValueServiceActor.{KvFailure, KvGet, KvResponse}

import scala.concurrent.Future
import scala.util.{Failure, Success}

abstract class KeyValueReadActor(override val threshold: Int, override val serviceRegistryActor: ActorRef)
    extends EnhancedThrottlerActor[CommandAndReplyTo[KvGet]] {
  override def receive = enhancedReceive.orElse(super.receive)

  override def processHead(head: CommandAndReplyTo[KvGet]) = instrumentedProcess {
    val processed = processGet(head.command)
    processed onComplete {
      case Success(response) => head.replyTo ! response
      case Failure(f) => head.replyTo ! KvFailure(head.command, f)
    }
    // This method should return how many "operations" have been performed to enable instrumentation of throughput
    // In this case it's 1 because we processed 1 KvGet
    processed.map(_ => 1)
  }

  def processGet(get: KvGet): Future[KvResponse]

  override protected lazy val instrumentationPath =
    KeyValueServiceActor.InstrumentationPath.concatNel(NonEmptyList.one("read"))
  override protected lazy val instrumentationPrefix = InstrumentationPrefixes.ServicesPrefix
  override def commandToData(snd: ActorRef) = { case get: KvGet =>
    CommandAndReplyTo(get, snd)
  }
}
