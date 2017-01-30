package cromwell.services.keyvalue

import akka.actor.{Actor, ActorRef}
import cromwell.core.{JobKey, WorkflowId}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{ExecutionContextExecutor, Future}
import scala.util.{Failure, Success}

object KeyValueServiceActor {

  final case class KvJobKey(callFqn: String, callIndex: Option[Int], callAttempt: Int)
  object KvJobKey {
    def apply(jobKey: JobKey): KvJobKey = KvJobKey(jobKey.scope.fullyQualifiedName, jobKey.index, jobKey.attempt)
  }

  final case class ScopedKey(workflowId: WorkflowId, jobKey: KvJobKey, key: String)

  sealed trait KvMessage {
    def key: ScopedKey
  }
  sealed trait KvMessageWithAction extends KvMessage {
    val action: KvAction
    def key = action.key
  }

  sealed trait KvAction extends KvMessage with ServiceRegistryMessage { override val serviceName = "KeyValue" }

  final case class KvPut(pair: KvPair) extends KvAction { override def key = pair.key }
  final case class KvGet(key: ScopedKey) extends KvAction

  sealed trait KvResponse extends KvMessage

  final case class KvPair(key: ScopedKey, value: Option[String]) extends KvResponse
  final case class KvFailure(action: KvAction, failure: Throwable) extends KvResponse with KvMessageWithAction
  final case class KvKeyLookupFailed(action: KvGet) extends KvResponse with KvMessageWithAction
  final case class KvPutSuccess(action: KvPut) extends KvResponse with KvMessageWithAction
}

trait KeyValueServiceActor extends Actor {
  implicit val ec: ExecutionContextExecutor

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
