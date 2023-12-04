package cromwell.services.keyvalue

import akka.actor.SupervisorStrategy.Escalate
import akka.actor.{Actor, ActorInitializationException, ActorLogging, OneForOneStrategy, Props}
import cats.data.NonEmptyList
import cromwell.core.{JobKey, WorkflowId}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

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

  final case class KvPair(key: ScopedKey, value: String) extends KvResponse
  final case class KvFailure(action: KvAction, failure: Throwable) extends KvResponse with KvMessageWithAction
  final case class KvKeyLookupFailed(action: KvGet) extends KvResponse with KvMessageWithAction
  final case class KvPutSuccess(action: KvPut) extends KvResponse with KvMessageWithAction

  val InstrumentationPath = NonEmptyList.of("keyvalue")
}

trait KeyValueServiceActor extends Actor with GracefulShutdownHelper with ActorLogging {
  protected def kvReadActorProps: Props
  protected def kvWriteActorProps: Props

  override val supervisorStrategy = OneForOneStrategy() {
    case _: ActorInitializationException => Escalate
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }

  private val kvReadActor = context.actorOf(kvReadActorProps, "KvReadActor")
  private val kvWriteActor = context.actorOf(kvWriteActorProps, "KvWriteActor")

  override def receive = {
    case get: KvGet => kvReadActor forward get
    case put: KvPut => kvWriteActor forward put
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.one(kvWriteActor))
  }
}
