package cromwell.services.keyvalue.impl

import akka.actor.{Props, ActorRef, Actor}
import com.typesafe.config.Config
import cromwell.database.CromwellDatabase
import cromwell.services.keyvalue.KeyValueService._

import scala.concurrent.Future
import scala.util.{Failure, Success}

object KeyValueServiceActor {
  def props(serviceConfig: Config, globalConfig: Config) = {
    Props(KeyValueServiceActor(serviceConfig, globalConfig))
  }
}

case class KeyValueServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with CromwellDatabase {
  private implicit val ec = context.dispatcher
  private implicit val system = context.system
  private var store: Map[ScopedKey, Option[String]] = Map.empty

  def receive = {
    case action: KvGet => respond(sender(), action, doGet(action))
    case action: KvPut => respond(sender(), action, doPut(action))
  }

  private def respond(replyTo: ActorRef, action: KvAction, response: Future[KvResponse]): Unit = {
    response.onComplete {
      case Success(x) => replyTo ! x
      case Failure(ex) => replyTo ! KvFailure(action, ex)
    }
  }

  private def doPut(put: KvPut): Future[KvResponse] = Future.successful {
    store += put.pair.key -> put.pair.value
    KvPutSuccess(put)
  }

  private def doGet(get: KvGet): Future[KvResponse] = Future.successful {
    store.get(get.key) match {
      case Some(maybeValue) => KvPair(get.key, maybeValue)
      case None => KvKeyLookupFailed(get)
    }
  }
}
