package cromwell.backend

import akka.actor.{Actor, ActorLogging}
import cromwell.backend.BackendActor._
import cromwell.backend.model.Subscription

object BackendActor {
  sealed trait BackendActorMessage
  case object Prepare extends BackendActorMessage
  case object Execute extends BackendActorMessage
  case object Stop extends BackendActorMessage
  case object CleanUp extends BackendActorMessage
  case class SubscribeToEvent[T](subscription: Subscription[T]) extends BackendActorMessage
  case class UnsubscribeToEvent[T](subscription: Subscription[T]) extends BackendActorMessage
}

/**
  * Defines basic structure and functionality to interact with a
  * Backend through an Akka actor system.
  * Backend methods should be implemented by each custom backend.
  */
trait BackendActor extends Backend with Actor with ActorLogging {
  def receive: Receive = {
    case Prepare => prepare
    case Execute => execute
    case Stop => stop
    case CleanUp => cleanUp
    case SubscribeToEvent(obj) => subscribeToEvent(obj)
    case UnsubscribeToEvent(obj) => unsubscribeToEvent(obj)
  }
}