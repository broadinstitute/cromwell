package cromwell.services.keyvalue

import akka.actor.{Actor, ActorRef}
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{ExecutionContextExecutor, Future}
import scala.util.{Failure, Success}


object KeyValueServiceActor {
  sealed  trait KvMessage

  sealed trait KvAction extends KvMessage with ServiceRegistryMessage {
    def serviceName = "KeyValue"
  }
  case class KvJobKey(callFqn: String, callIndex: Option[Int], callAttempt: Int)
  case class ScopedKey(workflowId: WorkflowId, jobKey: KvJobKey, key: String)
  case class KvPut(pair: KvPair) extends KvAction
  case class KvGet(key: ScopedKey) extends KvAction

  sealed trait KvResponse extends KvMessage
  case class KvPair(key: ScopedKey, value: Option[String]) extends KvResponse
  case class KvFailure(action: KvAction, failure: Throwable) extends KvResponse
  case class KvKeyLookupFailed(action: KvGet) extends KvResponse
  case class KvPutSuccess(action: KvPut) extends KvResponse
}

trait KeyValueServiceActor extends Actor {
  implicit val ec: ExecutionContextExecutor
  val serviceConfig: Config
  val globalConfig: Config

  def receive = {
    case action: KvGet => respond(sender(), action, doGet(action))
    case action: KvPut => respond(sender(), action, doPut(action))
  }

  def doPut(put: KvPut): Future[KvResponse]

  def doGet(get: KvGet): Future[KvResponse]

  private def respond(replyTo: ActorRef, action: KvAction, response: Future[KvResponse]): Unit = {
    response.onComplete {
      case Success(x) => replyTo ! x
      case Failure(ex) => replyTo ! KvFailure(action, ex)
    }
  }
}
