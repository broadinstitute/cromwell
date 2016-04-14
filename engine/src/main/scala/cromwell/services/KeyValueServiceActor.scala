package cromwell.services

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend.BackendJobDescriptor
import cromwell.engine.db.DataAccess
import cromwell.services.KeyValueServiceActor._
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

import scala.concurrent.{Await, Future}
import scala.util.{Failure, Success}

object KeyValueServiceActor {
  trait KvMessage
  trait KvAction extends KvMessage with ServiceRegistryMessage {
    def serviceName = "KeyValue"
  }
  trait KvResponse extends KvMessage

  case class ScopedKey(jobDescriptor: BackendJobDescriptor, key: String)
  case class KvPut(pair: KvPair) extends KvAction
  case class KvGet(key: ScopedKey) extends KvAction

  case class KvPair(key: ScopedKey, value: Option[String]) extends KvResponse
  case class KvFailure(action: KvAction, failure: Throwable) extends KvResponse
  case class KvKeyLookupFailed(action: KvGet) extends KvResponse
  case class KvPutSuccess(action: KvPut) extends KvResponse

  def props(serviceConfig: Config, globalConfig: Config) = {
    Props(KeyValueServiceActor(serviceConfig, globalConfig))
  }
}

case class KeyValueServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor {
  private implicit val ec = context.dispatcher
  private val dataAccess = DataAccess.globalDataAccess

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
    dataAccess.upsertExecutionInfo(
      put.pair.key.jobDescriptor.descriptor.id,
      put.pair.key.jobDescriptor.key,
      Map(put.pair.key.key -> put.pair.value),
      context.system
    ).map(_ => KvPutSuccess(put))
  }

  private def doGet(get: KvGet): Future[KvResponse] = {
    val executionInfo = dataAccess.getExecutionInfoByKey(
      get.key.jobDescriptor.descriptor.id,
      get.key.jobDescriptor.key.call,
      get.key.jobDescriptor.key.attempt,
      get.key.key
    )

    executionInfo map {
      case Some(maybeValue) => KvPair(ScopedKey(get.key.jobDescriptor, get.key.key), maybeValue)
      case None => KvKeyLookupFailed(get)
    }
  }
}