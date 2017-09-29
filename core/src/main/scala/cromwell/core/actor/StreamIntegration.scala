package cromwell.core.actor

import akka.actor.ActorRef
import akka.stream.QueueOfferResult

object StreamIntegration {
  trait StreamContext {
    def replyTo: ActorRef
    def request: Any
    /*
      clientContext way of attaching some information to the request which will get handed back at the end. In
        a sense, this could be viewed as the moral equivalent of the client keeping its own Map tracking requests to
        this information.
      */
    def clientContext: Option[Any] = None
  }
  case class EnqueueResponse(response: QueueOfferResult, request: StreamContext)
  case class BackPressure(request: Any)
  case class FailedToEnqueue(failure: Throwable, request: StreamContext)
}
