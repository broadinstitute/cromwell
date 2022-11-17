package cromwell.services

import akka.actor.{Actor, ActorRef}
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.scalalogging.StrictLogging
import cromwell.services.ServiceRegistryActor.{IoActorRef, NoIoActorRefAvailable, RequestIoActorRef}

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.concurrent.duration._
import scala.util.{Failure, Success}

trait IoActorRequester extends StrictLogging { this: Actor =>

  val serviceRegistryActor: ActorRef

  var _ioActorPromise: Option[Promise[ActorRef]] = None

  def requestIoActor(backoffInterval: FiniteDuration = 1.minute): Future[ActorRef] = _ioActorPromise match {
    case Some(promise) => promise.future
    case None =>
      val newPromise = Promise[ActorRef]()
      _ioActorPromise = Option(newPromise)
      requestIoActorInner(newPromise, backoffInterval)
      newPromise.future
  }

  private def requestIoActorInner(promise: Promise[ActorRef], backoffInterval: FiniteDuration): Unit = {
    implicit val ec: ExecutionContext = context.system.dispatcher
    implicit val timeout: Timeout = new Timeout(1.minute)

    serviceRegistryActor ? RequestIoActorRef onComplete {
      case Success(IoActorRef(actorRef)) =>
        promise.complete(Success(actorRef))
      case Success(NoIoActorRefAvailable) =>
        logger.warn(s"No IoActorRef available for ${self.path} yet. Retrying in $backoffInterval.")
        context.system.scheduler.scheduleOnce(backoffInterval) { requestIoActorInner(promise, backoffInterval) }
      case Success(other) =>
        val message = s"Programmer Error: Unexpected response to a RequestIoActor message in ${self.path}'s IoActorRequester: $other"
        logger.error(message)
        promise.failure(new Exception(message))
      case Failure(reason) =>
        promise.failure(new Exception("Failed to get an IoActor reference", reason))
    }
  }
}
