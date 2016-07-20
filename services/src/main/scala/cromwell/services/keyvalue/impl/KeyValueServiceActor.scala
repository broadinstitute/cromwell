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

case class KeyValueServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with KeyValueDatabaseAccess with CromwellDatabase {
  private implicit val ec = context.dispatcher

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
    upsertExecutionInfo(
      put.pair.key.workflowId,
      put.pair.key.jobKey,
      Map(put.pair.key.key -> put.pair.value),
      context.system
    ).map(_ => KvPutSuccess(put))
  }

  private def doGet(get: KvGet): Future[KvResponse] = {
    val executionInfo = getExecutionInfoByKey(
      get.key.workflowId,
      get.key.jobKey.scope,
      get.key.jobKey.attempt,
      get.key.key
    )

    executionInfo map {
      case Some(maybeValue) => KvPair(get.key, maybeValue)
      case None => KvKeyLookupFailed(get)
    }
  }
}