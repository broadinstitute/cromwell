package cromwell.services.keyvalue

import akka.actor.ActorRef
import cromwell.core.actor.BatchActor.CommandAndReplyTo
import cromwell.core.actor.ThrottlerActor
import cromwell.core.{JobKey, WorkflowId}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.Future
import scala.util.{Failure, Success}

object KeyValueServiceActor {
  final case class KvJobKey(callFqn: String, callIndex: Option[Int], callAttempt: Int)
  object KvJobKey {
    def apply(jobKey: JobKey): KvJobKey = KvJobKey(jobKey.node.fullyQualifiedName, jobKey.index, jobKey.attempt)
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

trait KeyValueServiceActor extends ThrottlerActor[CommandAndReplyTo[KvAction]] {
  override def commandToData(snd: ActorRef): PartialFunction[Any, CommandAndReplyTo[KvAction]] = {
    case c: KvAction => CommandAndReplyTo(c, snd)
  }

  override def processHead(action: CommandAndReplyTo[KvAction]) = action.command match {
    case get: KvGet => respond(action.replyTo, get, doGet(get))
    case put: KvPut => respond(action.replyTo, put, doPut(put))
  }

  def doPut(put: KvPut): Future[KvResponse]
  def doGet(get: KvGet): Future[KvResponse]

  private def respond(replyTo: ActorRef, action: KvAction, response: Future[KvResponse]): Future[KvResponse] = {
    response.onComplete {
      case Success(x) => replyTo ! x
      case Failure(ex) => replyTo ! KvFailure(action, ex)
    }
    response
  }
}
