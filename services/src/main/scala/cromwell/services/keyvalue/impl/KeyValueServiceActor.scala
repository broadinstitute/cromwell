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

case class KeyValueServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with BackendKeyValueDatabaseAccess with CromwellDatabase {
  private implicit val ec = context.dispatcher
  private implicit val system = context.system
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

  private def doPut(put: KvPut): Future[KvResponse] = {
    put.pair.value match {
      case Some(backendVal) => updateBackendKeyValuePair(put.pair.key.workflowId,
        put.pair.key.jobKey,
        put.pair.key.key,
        put.pair.value.get).map(_ => KvPutSuccess(put))
      case None => Future(KvFailure(put, new RuntimeException(s"Failed to find the value associated to key: ${put.pair.key.key}. This key cannot be added to the BackendKVStore.")))
    }
  }

  private def doGet(get: KvGet): Future[KvResponse] = {
    val backendValue = getBackendValueByKey(
      get.key.workflowId,
      get.key.jobKey.scope,
      get.key.jobKey.index,
      get.key.jobKey.attempt,
      get.key.key
    )

    backendValue map {
      case Some(maybeValue) => KvPair(get.key, Option(maybeValue))
      case None => KvKeyLookupFailed(get)
    }
  }
}
