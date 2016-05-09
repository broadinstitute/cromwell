package cromwell.util

import akka.actor._

import scala.concurrent.{Future, Promise, ExecutionContext}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object PromiseActor {
  /**
    * Sends a message to an actor and returns the future associated with the fullfilment of the reply
    * Can be used instead of the akka `ask` semantics, without any timeout
    *
    * @param sendTo      The actor that is wished to be asked something
    * @param message     The message to send to the actor
    * @param actorSystem ActorSystem from which to create the actor from
    * @return The future associated with the answer to the `ask`
    */
  private def askNoTimeout(message: Any, sendTo: ActorRef)(implicit actorSystem: ActorSystem): Future[Any] = {
    val promise = Promise[Any]()
    val promiseActor = actorSystem.actorOf(props(promise))
    sendTo.tell(message, promiseActor)
    promise.future
  }

  def props(promise: Promise[Any]): Props = Props(new PromiseActor(promise))

  implicit class EnhancedActorRef(val actorRef: ActorRef) extends AnyVal {
    def askNoTimeout(message: Any)(implicit ec: ExecutionContext, actorSystem: ActorSystem): Future[Any] = {
      PromiseActor.askNoTimeout(message, actorRef)
    }
  }
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
private class PromiseActor(promise: Promise[Any]) extends Actor {
  override def receive = {
    case Status.Failure(f) =>
      promise.tryFailure(f)
      context.stop(self)
    case success =>
      promise.trySuccess(success)
      context.stop(self)
  }
}
