package cromwell.services.keyvalue

import akka.actor.{ActorRef, Props}
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

final class InMemoryKvServiceActor extends KeyValueServiceActor {
  implicit val ec: ExecutionContext = context.dispatcher

  var kvStore = Map.empty[ScopedKey, String]

  override def receive = {
    case get: KvGet => respond(sender(), get, doGet(get))
    case put: KvPut => respond(sender(), put, doPut(put))
  }

  def doGet(get: KvGet): Future[KvResponse] = kvStore.get(get.key) match {
    case Some(value) => Future.successful(KvPair(get.key, value))
    case None => Future.successful(KvKeyLookupFailed(get))
  }

  def doPut(put: KvPut): Future[KvResponse] = {
    kvStore += (put.key -> put.pair.value)
    Future.successful(KvPutSuccess(put))
  }

  override protected def kvReadActorProps = Props.empty

  override protected def kvWriteActorProps = Props.empty

  private def respond(replyTo: ActorRef, action: KvAction, response: Future[KvResponse]): Unit = {
    response.onComplete {
      case Success(x) => replyTo ! x
      case Failure(ex) => replyTo ! KvFailure(action, ex)
    }
  }
}
