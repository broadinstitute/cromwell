package cromwell.util

import akka.actor._

import scala.concurrent.{Future, Promise}

private class PromiseActor(promise: Promise[Any], sendTo: ActorRef, msg: Any) extends Actor with ActorLogging {

  context.watch(sendTo)
  sendTo ! msg

  override def receive = {
    case Status.Failure(f) =>
      promise.tryFailure(f)
      context.stop(self)
    case Terminated(actorRef) =>
      if (actorRef == sendTo) {
        promise.tryFailure(new RuntimeException("Promise-watched actor completed before sending back a message"))
      } else {
        log.error("Spooky happenstances! A Terminated({}) message  was sent to a private Promise actor which wasn't watching it!?", actorRef)
      }
      context.stop(self)
    case success =>
      promise.trySuccess(success)
      context.stop(self)
  }
}

object PromiseActor {
  /**
    * Sends a message to an actor and returns the future associated with the fullfilment of the reply
    * Can be used instead of the akka `ask` semantics, without any timeout
    *
    * @param sendTo          The actor that is wished to be asked something
    * @param message         The message to send to the actor
    * @param actorRefFactory Factory from which to create the actor from
    * @return The future associated with the answer to the `ask`
    */
  private def askNoTimeout(message: Any, sendTo: ActorRef)(implicit actorRefFactory: ActorRefFactory): Future[Any] = {
    val promise = Promise[Any]()
    val _ = actorRefFactory.actorOf(props(promise, sendTo, message))
    promise.future
  }

  def props(promise: Promise[Any], sendTo: ActorRef, msg: Any): Props = Props(new PromiseActor(promise, sendTo, msg))

  implicit class EnhancedActorRef(val actorRef: ActorRef) extends AnyVal {
    def askNoTimeout(message: Any)(implicit actorRefFactory: ActorRefFactory): Future[Any] = {
      PromiseActor.askNoTimeout(message, actorRef)
    }
  }
}
