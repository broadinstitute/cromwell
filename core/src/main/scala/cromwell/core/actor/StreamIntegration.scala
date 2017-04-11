package cromwell.core.actor

import akka.actor.ActorRef
import akka.stream.QueueOfferResult

object StreamIntegration {
  trait StreamContext {
    def replyTo: ActorRef
    def request: Any
    def clientContext: Option[Any] = None
  }
  case class EnqueueResponse(response: QueueOfferResult, request: StreamContext)
  case class BackPressure(request: Any)
  case class FailedToEnqueue(failure: Throwable, request: StreamContext)
}
