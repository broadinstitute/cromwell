package cromwell.services.pubsub

import akka.actor.Actor
import com.typesafe.scalalogging.LazyLogging
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.pubsub.PubSubServiceActor.Publish
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

trait PubSubServiceActor extends Actor with LazyLogging {
  override def receive: Receive = {
    case Publish(topic, message) =>
      logger.debug("Pubsub service received request to publish message " + message + " onto topic " + topic)
      publish(topic, message)
    case ShutdownCommand => context.stop(self) // Not necessary but service registry requires it. See #2575
  }

  protected def publish(topic: String, message: String): Unit
}

object PubSubServiceActor {
  sealed abstract class PubSubServiceActorMessage extends ServiceRegistryMessage { override def serviceName = "PubSub" }
  case class Publish(topic: String, message: String) extends PubSubServiceActorMessage
}